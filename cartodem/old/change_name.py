import sqlite3
fname = 'Bhuvan_CartoDEM_v3r1_TerrainRGB.mbtiles'
name = fname.replace('.mbtiles', '')
con = sqlite3.connect(fname)
#res = con.execute('SELECT * FROM sqlite_schema')
#for r in res:
#    print(r)
#print(res)
{'attribution': '', 'description': '', 'format': 'png', 'name': 'Countries'}
kv = {
    'name': name,
    'maxzoom': '12',
    'minzoom': '5',
    'attribution': '<a href="https://bhuvan-app3.nrsc.gov.in/data/download/index.php" target="_blank" rel="noopener noreferrer">Bhuvan CartoDEM 30m v3r1</a>',
    'format': 'webp',
    'type': 'baselayer',
    'description': '',
    'version': '1',

}
#sql = f"DELETE from metadata WHERE name='{name}'"
#con.execute(sql)
#con.commit()
#exit(0)
for k, v in kv.items():
    sql = f"INSERT into metadata (name, value) VALUES ('{k}','{v}');"
    #sql = f"UPDATE metadata SET value='{v}' WHERE name='{k}'"
    print(f'executing {sql}')
    con.execute(sql)
    con.commit()
