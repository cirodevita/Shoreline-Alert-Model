import sys
import json

#excludeWavemetersById = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 313, 315, 328, 455]
excludeWavemetersById = [69, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 313, 314, 315, 316, 328, 329, 455, 456]
features = []

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python " + str(sys.argv[0]) + " wavemeter.geojson")
        sys.exit(-1)

    wavemeterPath = sys.argv[1]
    with open(wavemeterPath, "r") as jsonFile:
        wavemeters = json.load(jsonFile)

    for idx, wavemeter in enumerate(wavemeters["features"]):
        id = wavemeter["properties"]["id"]
        if id not in excludeWavemetersById:
            features.append(wavemeter)

    wavemeters["features"] = features

    with open(wavemeterPath, "w") as jsonFile:
        json.dump(wavemeters, jsonFile)
