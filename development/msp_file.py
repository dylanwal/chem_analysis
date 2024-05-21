
import numpy as np

import chem_analysis as ca
import chem_analysis.utils.math as utils_math


class ChemicalCompound:
    def __init__(self, name, formula, mw, casno, id_, ri, comment, num_peaks, peaks):
        self.name = name
        self.formula = formula
        self.mw = mw
        self.casno = casno
        self.id = id_
        self.ri = ri
        self.comment = comment
        self.num_peaks = num_peaks
        self.peaks = peaks

    def __repr__(self):
        return (f"ChemicalCompound(name={self.name}, formula={self.formula}, mw={self.mw}, "
                f"casno={self.casno}, id={self.id}, ri={self.ri}, comment={self.comment}, "
                f"num_peaks={self.num_peaks}, peaks={self.peaks})")


def parse_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = {}
    i = 0
    peaks = None
    compounds = []
    for line in lines:
        if line == "\n":
            compound = ChemicalCompound(
                name=data["name"],
                formula=data["formula"],
                mw=data["mw"],
                casno=data["casno"],
                id_=data["id"],
                ri=data["ri"] if "ri" in data else None,
                comment=data["comment"],
                num_peaks=data["num_peaks"],
                peaks=peaks.astype(utils_math.min_int_dtype_array(peaks))
            )
            compounds.append(compound)
            data = {}
            peaks = None
            continue

        line = line.strip()
        if line.startswith("Name:"):
            data["name"] = line.split(":")[1].strip()
        elif line.startswith("Formula:"):
            data["formula"] = line.split(":")[1].strip()
        elif line.startswith("MW:"):
            data["mw"] = float(line.split(":")[1])
        elif line.startswith("CASNO:"):
            data["casno"] = line.split(":")[1].strip()
        elif line.startswith("ID:"):
            data["id"] = int(line.split(":")[1])
        elif line.startswith("RI:"):
            data["ri"] = int(line.split(":")[1])
        elif line.startswith("Comment:"):
            data["comment"] = line.split(":")[1].strip()
        elif line.startswith("Num peaks:"):
            data["num_peaks"] = int(line.split(":")[1])
            if peaks is not None:
                raise ValueError("peaks should be None.")
            peaks = np.empty((data["num_peaks"], 2), dtype="uint64")
            i = 0
        else:
            line = line.split(" ")
            if len(line) != 2:
                raise ValueError("Incorrect number of line")
            mz = int(line[0])
            intensity = int(float(line[1])*100)
            peaks[i, :] = (mz, intensity)
            i += 1

    if peaks:
        compound = ChemicalCompound(
            name=data["name"],
            formula=data["formula"],
            mw=data["mw"],
            casno=data["casno"],
            id_=data["id"],
            ri=data["ri"],
            comment=data["comment"],
            num_peaks=data["num_peaks"],
            peaks=peaks
        )
        compounds.append(compound)

    return compounds


def main():
    i = 4
    file_path = fr"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\NIST_20\NIST20_{i}.MSP"
    compounds = parse_file(file_path)

    lib = ca.gc_lc.GCLibrary(f"NIST20_{i}")
    for compound in compounds:
        comp = ca.gc_lc.Compound(
            name=compound.name + "|" + compound.formula,
            label=compound.name,
            cas=compound.casno,
            responses=[
                ca.gc_lc.CompoundResponse(
                    method="ms",
                    mass_spectrum=compound.peaks
                )
            ],
        )
        lib.add_compound(comp)

    lib.to_JSON(fr"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\nist20_{i}.json", overwrite=True, binary=True, optimize=True)


if __name__ == "__main__":
    main()
