"""Install checksum-pinned standalone WhiteboxTools binaries inside --workdir.

Uses the conda-forge Linux ARM64 build without creating a Conda environment or
installing its Python wrapper. A private libgcc runtime is supplied to the CLI.
Requires Python 3.14's standard-library Zstandard reader.
"""
import argparse
import compression.zstd
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tarfile
import urllib.request
import zipfile

PACKAGES = {
    'whitebox_tools-2.4.0-py314hfe60d44_2.conda': 'e3e680dc957a87fa99ed73f9cd0fc9028866557684af28b09ac9537e76dc2102',
    'libgcc-14.3.0-he277a41_7.conda': '2e3ed044b4c33e7d2633b13e25e3764053e1a4eca64886bb88a0de96a0e03dea',
}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    args=parser.parse_args();root=args.workdir.expanduser().resolve()
    if (platform.system(),platform.machine())!=('Linux','aarch64'):
        parser.error('These pinned packages target Linux ARM64')
    folder=root/'tools/whitebox-prebuilt';folder.mkdir(parents=True,exist_ok=True)
    records=[]
    for name,expected in PACKAGES.items():
        package=folder/name
        url='https://conda.anaconda.org/conda-forge/linux-aarch64/'+name
        if not package.exists():
            temporary=package.with_suffix('.part')
            with urllib.request.urlopen(url) as response,temporary.open('wb') as out:
                shutil.copyfileobj(response,out)
            if digest(temporary)!=expected:raise ValueError('Package checksum mismatch: '+name)
            temporary.replace(package)
        if digest(package)!=expected:raise ValueError('Package checksum mismatch: '+name)
        with zipfile.ZipFile(package) as z:
            for part in z.namelist():
                if not part.endswith('.tar.zst'):continue
                # Keep both packages' license/recipe metadata without collisions.
                destination=folder if part.startswith('pkg-') else folder/'package-info'/name
                destination.mkdir(parents=True,exist_ok=True)
                with compression.zstd.ZstdFile(z.open(part)) as stream,tarfile.open(fileobj=stream,mode='r|') as tar:
                    tar.extractall(destination,filter='data')
        records.append(dict(url=url,sha256=expected))
    binary=folder/'bin/whitebox_tools'
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(folder/'lib')+(':'+env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
    settings=dict(working_directory='',verbose_mode=True,max_procs=1,compress_rasters=True)
    (binary.parent/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    version=subprocess.check_output([str(binary),'--version'],env=env,text=True).strip()
    for tool in ('CostDistance','CostAllocation','BreachDepressions','D8Pointer','D8FlowAccumulation'):
        subprocess.run([str(binary),'--toolhelp='+tool],env=env,check=True,stdout=subprocess.DEVNULL)
    runtime=folder/'lib/libgcc_s.so.1'
    record=dict(distribution='conda-forge standalone CLI',version=version,packages=records,
                binary=str(binary.relative_to(root)),binary_sha256=digest(binary),
                library_directory=str((folder/'lib').relative_to(root)),
                runtime_sha256={'libgcc_s.so.1':digest(runtime)},settings=settings,
                installer_sha256=digest(Path(__file__)),
                license='WhiteboxTools: MIT; libgcc: GPL-3.0-only WITH GCC-exception-3.1',
                compiled_locally=False)
    (folder/'build.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__=='__main__':main()
