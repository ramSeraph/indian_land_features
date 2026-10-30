"""Build the pinned WhiteboxTools open core inside a chosen workspace.

No system installation or shell-profile changes are made. Rust, Cargo caches,
sources, build products and provenance all stay beneath --workdir/tools.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import tarfile
import time
import urllib.request

REVISION = '28ed01f506f2024cfb52db85348b887f2ded4e5d'
SOURCE_SHA256 = '898e24b73c5a5145071029c558f465f58b1742d0e15f55f77cd5fe112c347951'
LOCK_SHA256 = '582980e7fe4504499ae95e991ed48d46927034510c694da7be61e69750834c41'
CORE_LOCK_SHA256 = 'ff22aa0cf7aae0eb65bd2c80d1261ab5519b4c47a8e16c41fd8a3e71b6180131'
RUST_VERSION = '1.79.0'  # Compatible with upstream's locked time 0.3.20.
RUSTUP_VERSION = '1.28.2'
RUSTUP_SHA256 = {
    'aarch64-unknown-linux-gnu': 'e3853c5a252fca15252d07cb23a1bdd9377a8c6f3efa01531109281ae47f841c',
}
HYDRO_TOOLS = {
    'FlowAccumulationFullWorkflow': 'flow_accum_full_workflow',
    'Watershed': 'watershed', 'Basins': 'basins',
    'D8Pointer': 'd8_pointer', 'D8FlowAccumulation': 'd8_flow_accum',
    'BreachDepressions': 'breach_depressions',
    'BreachDepressionsLeastCost': 'breach_depressions_least_cost',
    'FillDepressionsWangAndLiu': 'fill_depressions_wang_and_liu',
}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def download(url, target, expected):
    if not target.exists():
        temporary = target.with_suffix(target.suffix + '.part')
        with urllib.request.urlopen(url) as response, temporary.open('wb') as stream:
            shutil.copyfileobj(response, stream)
        if digest(temporary) != expected:
            raise ValueError(f'Checksum mismatch for {target.name}')
        temporary.replace(target)
    if digest(target) != expected:
        raise ValueError(f'Checksum mismatch for {target.name}')


def focus_hydrology(archive, source):
    """Change only tool registration; verify every selected algorithm unchanged."""
    manager = 'whitebox-tools-app/src/tools/mod.rs'
    hydro = 'whitebox-tools-app/src/tools/hydro_analysis/mod.rs'
    hashes = {}
    with tarfile.open(archive) as tar:
        def upstream(name):
            return tar.extractfile(f'whitebox-tools-{REVISION}/'+name).read()
        text = upstream(manager).decode()
        text = re.sub(r'^pub mod \w+;\n', '', text, flags=re.M)
        text = 'pub mod hydro_analysis;\n'+text
        start = text.index('        let mut tool_names = vec![];')
        end = text.index('        tool_names.sort();', start)
        names = '        let mut tool_names = vec![];\n'+''.join(
            f'        tool_names.push("{name}".to_string());\n' for name in HYDRO_TOOLS)
        text = text[:start]+names+text[end:]
        start = text.index('    fn get_tool(&self,')
        end = text.index('    fn get_plugin_list(', start)
        method = "    fn get_tool(&self, tool_name: &str) -> Option<Box<dyn WhiteboxTool + 'static>> {\n"
        method += '        match tool_name.to_lowercase().replace("_", "").as_ref() {\n'
        method += ''.join(f'            "{name.lower()}" => Some(Box::new(hydro_analysis::{name}::new())),\n'
                          for name in HYDRO_TOOLS)
        method += '            _ => None,\n        }\n    }\n\n'
        (source/manager).write_text(text[:start]+method+text[end:])
        (source/hydro).write_text(''.join(f'mod {module};\npub use self::{module}::{name};\n'
                                        for name,module in HYDRO_TOOLS.items()))
        for module in HYDRO_TOOLS.values():
            name=f'whitebox-tools-app/src/tools/hydro_analysis/{module}.rs'
            expected=hashlib.sha256(upstream(name)).hexdigest()
            if digest(source/name)!=expected:
                raise ValueError(f'Upstream algorithm changed: {module}')
            hashes[module]=expected
    return dict(algorithms_sha256=hashes, registry_sha256={name:digest(source/name) for name in (manager,hydro)})


def compile_bounded(command, source, env, max_rss_mib=1400):
    """Limit aggregate compiler-process memory, including Cargo's descendants."""
    process=subprocess.Popen(command,cwd=source,env=env,start_new_session=True)
    peak=0;started=time.monotonic()
    def rss_tree(pid):
        try:
            folder=Path('/proc')/str(pid)
            fields=dict(line.split(':',1) for line in (folder/'status').read_text().splitlines() if ':' in line)
            rss=int(fields.get('VmRSS','0 kB').split()[0])
            # Cargo starts compiler children from worker threads, not just its
            # main thread. Linux records children separately for each task.
            children=set()
            for task in (folder/'task').iterdir():
                try:children.update(int(child) for child in (task/'children').read_text().split())
                except FileNotFoundError:pass
            return rss+sum(rss_tree(child) for child in children)
        except FileNotFoundError:
            return 0
    try:
        while process.poll() is None:
            peak=max(peak,rss_tree(process.pid))
            if peak>max_rss_mib*1024:
                raise RuntimeError(f'Compiler exceeded {max_rss_mib} MiB RSS limit')
            time.sleep(.2)
        if process.returncode:
            raise subprocess.CalledProcessError(process.returncode,command)
    finally:
        if process.poll() is None:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait()
    return dict(peak_observed_rss_mib=round(peak/1024,2),elapsed_seconds=round(time.monotonic()-started,2),
                memory_limit_mib=max_rss_mib,scope='process tree including children of worker threads')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != 'Linux' or platform.machine() != 'aarch64':
        parser.error('This pinned bootstrap currently supports Linux aarch64; other builds need a verified Rust bootstrap checksum')
    target = 'aarch64-unknown-linux-gnu'
    root = args.workdir.expanduser().resolve()
    folder = root / 'tools' / 'whitebox-2.4.0'
    folder.mkdir(parents=True, exist_ok=True)
    archive = folder / 'source.tar.gz'
    download(f'https://codeload.github.com/jblindsay/whitebox-tools/tar.gz/{REVISION}', archive, SOURCE_SHA256)
    source = folder / f'whitebox-tools-{REVISION}'
    if not source.exists():
        with tarfile.open(archive) as tar:
            tar.extractall(folder, filter='data')
    if digest(source / 'Cargo.lock') not in (LOCK_SHA256, CORE_LOCK_SHA256):
        raise ValueError('Cargo.lock changed')
    # Only compile the command-line core. The excluded plugin workspace has a
    # stale tsp-rs checksum; do not replace or accept that unrelated dependency.
    # The focused lock retains the original 141 core dependency package versions.
    core_lock = Path(__file__).with_name('whitebox-core.Cargo.lock')
    if digest(core_lock) != CORE_LOCK_SHA256:
        raise ValueError('Focused core lockfile changed')
    manifest = source / 'Cargo.toml'
    manifest.write_text(manifest.read_text().replace('"whitebox-plugins", ', '').replace('"whitebox-runner", ', ''))
    shutil.copy2(core_lock, source / 'Cargo.lock')
    focused=focus_hydrology(archive,source)

    bootstrap = folder / 'rustup-init'
    download(f'https://static.rust-lang.org/rustup/archive/{RUSTUP_VERSION}/{target}/rustup-init',
             bootstrap, RUSTUP_SHA256[target])
    bootstrap.chmod(0o755)
    env = os.environ.copy()
    env.update(CARGO_HOME=str(folder / 'cargo'), RUSTUP_HOME=str(folder / 'rustup'),
               CARGO_BUILD_JOBS='1', CARGO_TARGET_DIR=str(folder / 'target'))
    cargo = folder / 'cargo' / 'bin' / 'cargo'
    if not cargo.exists():
        subprocess.run([str(bootstrap), '-y', '--no-modify-path', '--profile', 'minimal',
                        '--default-toolchain', RUST_VERSION], env=env, check=True)
    env['PATH'] = str(cargo.parent) + os.pathsep + env['PATH']
    compilation=compile_bounded([str(cargo), '+' + RUST_VERSION, 'build', '--release', '--locked',
                    '--jobs', '1', '--package', 'whitebox_tools'], source, env)
    binary = folder / 'target' / 'release' / 'whitebox_tools'
    # The upstream default can use all processors; the bounded pilot uses one.
    (binary.parent / 'settings.json').write_text(json.dumps(dict(
        working_directory='', verbose_mode=True, max_procs=1,
        compress_rasters=True), indent=2) + '\n')
    record = dict(source_revision=REVISION, source_sha256=SOURCE_SHA256,
                  upstream_cargo_lock_sha256=LOCK_SHA256,
                  cargo_lock_sha256=digest(source / 'Cargo.lock'), rust_version=RUST_VERSION,
                  workspace_scope='command-line core; plugins and GUI excluded; dependency versions unchanged',
                  registered_tools=list(HYDRO_TOOLS), source_verification=focused,
                  compilation=compilation, installer_sha256=digest(Path(__file__)),
                  rustup_version=RUSTUP_VERSION, target=target,
                  binary_sha256=digest(binary), binary=str(binary.relative_to(root)),
                  license='MIT', jobs=1,
                  version=subprocess.check_output([str(binary), '--version'], text=True).strip())
    (folder / 'build.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()
