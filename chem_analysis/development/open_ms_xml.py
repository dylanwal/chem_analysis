
import xml.etree.ElementTree as ET
import base64
import dataclasses

import numpy as np

MAX_MASS_RANGE = 600


@dataclasses.dataclass
class MZData:
    sample_name: str
    acquisition_method: str
    spectrums: np.ndarray
    times: np.ndarray


def decoding_mzdata(file_path):
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
        mz = np.frombuffer(base64.b64decode(spec.find("mzArrayBinary").find("data").text), dtype="float64").astype("uint32")
        intensity = np.frombuffer(base64.b64decode(spec.find("intenArrayBinary").find("data").text), dtype="float32")
        spectrums[i, mz] = intensity

    data["spectrums"] = spectrums
    data["times"] = times

    return MZData(**data)


def main():
    file_path = r"C:\Users\nicep\Desktop\research_wis\data\10\10_13\DJW-10-13-600min-TMS.mzdata.xml"
    data = decoding_mzdata(file_path)

    x = data.times
    y = data.spectrums.sum(axis=1)
    import plotly.graph_objs as go
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y))
    fig.show()
    print(data)


if __name__ == "__main__":
    main()
