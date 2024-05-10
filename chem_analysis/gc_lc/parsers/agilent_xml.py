import pathlib
import xml.etree.ElementTree as ET
import base64
from typing import Any

import numpy as np

MAX_MASS_RANGE = 600


def decoding_mzdata(file_path: str | pathlib.Path) -> dict[str, Any]:
    with open(file_path, 'r') as file:
        xml_string = file.read()

    root = ET.fromstring(xml_string)

    data = {
        "sample_name": root.find("description").find("admin").find("sampleName").text,
        "acquisition_method": root.find("description").find("instrument").find("additional").find("userParam").get("value")
    }

    spectrum_list = root.find('spectrumList')
    spectrums = np.zeros((int(spectrum_list.get('count')), MAX_MASS_RANGE), dtype="float32")
    times = np.empty(int(spectrum_list.get('count')), dtype="float32")

    for i, spec in enumerate(spectrum_list):
        times[i] = np.array([spec.find("spectrumDesc").find("spectrumSettings").find("spectrumInstrument")[2].get("value")], dtype="float32")
        mz = np.round(np.frombuffer(base64.b64decode(spec.find("mzArrayBinary").find("data").text), dtype="float64")).astype("uint32")
        intensity = np.frombuffer(base64.b64decode(spec.find("intenArrayBinary").find("data").text), dtype="float32")
        if check_for_duplicates(mz):
            print('duplicate mz detected')
        spectrums[i, mz] = intensity

    data["spectrums"] = spectrums
    data["times"] = times

    return data


def check_for_duplicates(arr):
    unique_elements, counts = np.unique(arr, return_counts=True)
    if np.any(counts > 1) >= 2:
        return True
    return False
