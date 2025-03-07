# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "requests",
#     "imgcat",
#     "pillow",
#     "bs4"
# ]
# ///

import json
import base64
import hashlib

from io import BytesIO
from pathlib import Path
from urllib.parse import urljoin

import requests
from imgcat import imgcat
from PIL import Image
from bs4 import BeautifulSoup

B64_PREFIX = 'data:image/png;base64,'
MAIN_URL ='https://bhuvan-app3.nrsc.gov.in/data/download/index.php'

raw_dir = Path('data/raw/v3_r1/')

def get_tile_names():
    tnames = []
    with open('data/cdem_v3_r1.geojsonl', 'r') as f:
        for line in f:
            feat = json.loads(line)
            props = feat['properties']
            tname = props['ref_1'].lower()
            tnames.append(tname)
    return tnames

if __name__ == '__main__':
    tnames = get_tile_names()
    session = requests.session()
    resp = session.get(MAIN_URL)
    if not resp.ok:
        raise Exception('unable to get main page')

    resp_text = resp.text
    #Path('temp.html').write_text(resp_text)

    login_url = urljoin(MAIN_URL, 'login.php')
    resp = session.get(login_url, headers = { 'Referer': MAIN_URL })
    if not resp.ok:
        raise Exception('unable to get login page')
    resp_text = resp.text
    soup = BeautifulSoup(resp_text, 'html.parser')
    inp = soup.find('input', { 'name': 'execution' })
    execution = inp.get('value')
    #Path('temp1.html').write_text(resp_text)
    curr_url = resp.url

    captcha_url = urljoin(curr_url, '/cas/captcha/generate')
    print(captcha_url)
    resp = session.get(captcha_url, headers = { 'Referer': curr_url })
    if not resp.ok:
        raise Exception('Unable to get captcha')
    cdata = resp.json()
    cid = cdata['captchaId']
    img_b64 = cdata['image']
    img_b64 = img_b64[len(B64_PREFIX):]
    img_bytes = base64.b64decode(img_b64)
    imgcat(Image.open(BytesIO(img_bytes)))
    val = input('Enter Captcha: ')
    val = val.strip()
    user_data = json.loads(Path('data/user.json').read_text())
    username = user_data['name']
    password = user_data['password'].encode('utf-8')
    password = f'{hashlib.sha256(password)}&&{hashlib.md5(password)}'
    form_data = {
        'username': username,
        'password': password,
        'customFields[0]': cid,
        'customFields[1]': val,
        'execution': execution,
        '_eventId': 'submit',
        'geolocation': '',
    }
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
        'Referer': curr_url,
    }
    resp = session.post(login_url, data=form_data, headers=headers)
    if not resp.ok:
        print(resp.text)
        raise Exception('login failed')

    raw_dir.mkdir(exist_ok=True, parents=True)
    total = len(tnames)
    count = 0
    for tname in tnames:
        count += 1
        out_file = raw_dir / f'{tname}.zip'
        if out_file.exists():
            continue
        print(f'downloading {tname} - {count}/{total}')
        download_url = f'https://bhuvan-app3.nrsc.gov.in/isroeodatadownloadutility/tiledownloadnew_cfr_new.php?f=cdn{tname}_v3r1.zip&se=CDEM&u={username}'
        resp = session.get(download_url, headers = {'referer': MAIN_URL })
        if not resp.ok:
            print(resp.text)
            raise Exception('unable to download tile')
        out_file.write_bytes(resp.content)


    

