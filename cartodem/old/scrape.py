# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "requests",
# ]
# ///

import json
from pathlib import Path
from pprint import pprint
from datetime import datetime, timedelta

import requests
from urllib.parse import urljoin

BASE_URL = "https://bhoonidhi.nrsc.gov.in"
session_info = {"jwt": None,"userId": None, "user_email": None, "username": None, "password": None, "sid": None, "scenes": []}

BBOX = {
    'tllat': 40.403039,
    'tllon': 53.839874,
    'brlat': 2.578170,
    'brlon': 103.761749,
}

def login(username, password, req_session):
    login_url = urljoin(BASE_URL, "bhoonidhi/LoginServlet")
    payload = {
        "userId": username,
        "password": password,
        "action": "VALIDATE_LOGIN",
        "oldDB": "false"
    }
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    response = req_session.post(login_url, data=json.dumps(payload), headers=headers)
    if response.status_code == 200:
        response_data = response.json()["Results"][0]
        if "JWT" in response_data and response_data["JWT"]:
            return response_data
        else:
            print("Error logging in. \nReason:", response_data['MSG'])
    else:
        print("Error logging in. Status code:", response.status_code)
         
def validate_session(jwt, req_session):
    login_url = urljoin(BASE_URL, "bhoonidhi/LoginServlet")
    payload = {
        "action": "VALIDATE_SESSION",
    }
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "token": jwt
    }
    response = req_session.post(login_url, data=json.dumps(payload), headers=headers)
    if response.status_code == 200:
        print("Session validation successful")
    else:
        print("Error logging in. Status code:", response.status_code)
 
def authenticate(session_file, req_session):
    user_data = json.loads(Path('data/user.json').read_text())
    username = user_data['name']
    password = user_data['password']
    session = login(username, password, req_session)
    session_info["jwt"] = session["JWT"]
    session_info["userId"] = session["USERID"]
    session_info["user_email"] = session["USEREMAIL"]
    session_info["user_name"] = session["USERNAME"]
    session_info["username"] = username
    session_info["password"] = password
    validate_session(session_info["jwt"], req_session)
    session_file.write_text(json.dumps(session_info))

def get_scenes(req_session):
    start_date = datetime.now() - timedelta(days=30)
    end_date = datetime.now()
    sdate = start_date.strftime("%b%%2F%d%%2F%Y").upper()
    edate = end_date.strftime("%b%%2F%d%%2F%Y").upper()

    product = "Standard"
    sat_sen = "CartoSat-1_PAN_CartoDEM-30m"


    payload = {
        "userId": session_info.get('userId'),
        "prod": product,
        "selSats": sat_sen,
        "offset": "0",
        "sdate": sdate,
        "edate": edate,
        "query": "area",
        "queryType": "polygon",
        "isMX": "No",
        "filters": "%7B%7D"
    }
    payload.update(BBOX)

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "token": session_info.get("jwt")
    }
    url = urljoin(BASE_URL, "bhoonidhi/ProductSearch")
    response = req_session.post(url, headers=headers, json=payload)
    if response.status_code == 200:
        data = response.json()
        results = data["Results"]
        return results
    else:
        print("Request failed. Status code:", response.status_code)

def download_scene(scene, req_session):
    product_id = scene['ID']
    zname = f"{product_id}.zip"
    out_file = Path('data/raw/bhoonidhi') / zname 
    if out_file.exists():
        return
    cdate = datetime.now().strftime("%d %B %Y").upper()
    url = urljoin(BASE_URL, f'bhoonidhi/data/CARTODEM/P5/PAN/30m/{product_id}.zip')
    params = {
        'token': session_info['jwt'],
        'product_id': product_id,
        'cartDate': cdate,
        'sid': scene['srt']
    }
    resp = req_session.get(url, params=params)
    if not resp.ok:
        print(resp.status_code)
        print(resp.text)
        raise Exception(f'Unable to get {product_id}')
    out_file.write_bytes(resp.content)

    #add_to_cart(scene, req_session)


if __name__ == '__main__':

    req_session = requests.session()
    Path('data').mkdir(exist_ok=True)
    session_file = Path('data/session.json')
    if not session_file.exists():
        authenticate(session_file, req_session)
    else:
        session_info = json.loads(session_file.read_text())
    
    scenes_file = Path('data/scenes.json')
    if scenes_file.exists():
        scenes = json.loads(scenes_file.read_text())
    else:
        scenes = get_scenes(req_session)
        scenes_file.write_text(json.dumps(scenes))

    lines = Path('done.txt').read_text().split()
    done = set([ l.strip() for l in lines if l.strip() != ''])
    count = 0
    total = len(scenes)
    for scene in scenes:
        count += 1
        if scene["ID"] in done:
            continue
        print(f'handling {scene["ID"]} {count}/{total}')
        download_scene(scene, req_session)





