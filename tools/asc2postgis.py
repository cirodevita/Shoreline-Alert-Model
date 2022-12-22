import sys
import os
import psycopg2


# Fix configuration
def getConfigFromAsc(filePath):
    file = open(filePath, 'r')
    lines = file.read().splitlines()[:6]
    cfg = {}

    for line in lines:
        element = line.split(" ")
        cfg[element[0]] = element[1]

    file.close()

    return cfg


def saveIntoDb(filePath, cfg, connection):
    x = float(list(cfg.items())[2][1])
    y = float(list(cfg.items())[3][1])
    step = float(list(cfg.items())[4][1])

    cur = connection.cursor()
    file = open(filePath, 'r')
    rows = file.read().splitlines()[6:]
    for idx_i, i in enumerate(rows[::-1]):
        columns = i.split(" ")
        columns = [ele for ele in columns if ele.strip()]
        for idx_j, j in enumerate(columns):
            X = x + (idx_j*step)
            Y = y + (idx_i*step)
            sql = "INSERT INTO public.points (geom) VALUES (ST_SetSRID(ST_MakePoint(%s, %s, %s), 4326))"
            # print(X, Y, j)
            cur.execute(sql, (X, Y, j))

    connection.commit()
    file.close()


if len(sys.argv) != 2:
    print("Usage: python " + str(sys.argv[0]) + " DTM.asc")
    sys.exit(-1)

#lidarPath = sys.argv[1]
dtmPath = sys.argv[1]

connection = psycopg2.connect(host='127.0.0.1', dbname='dtm', user='dtm', password='Dtm2022')

if os.path.isfile(dtmPath):
    print(dtmPath)
    
    cfg = getConfigFromAsc(dtmPath)    
    saveIntoDb(dtmPath, cfg, connection)

"""
if os.path.isdir(lidarPath):
    for file in os.listdir(lidarPath):
        if file.endswith("asc"):
            path = os.path.join(lidarPath, file)
            print(path)

            cfg = getConfigFromAsc(path)
            saveIntoDb(path, cfg, connection)
else:
    sys.exit("Path is not a directory!")
"""

connection.close()

