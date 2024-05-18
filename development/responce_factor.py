import glob

import chem_analysis as ca

# lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\new_library.JSON"
# chemistry_lib = ca.mass_spec.GCLibrary.from_JSON(lib_path)

REFERENCE_DATA = [
     {
            "label": "K6_2",
            "group": "ketone",
            "name": "2-hexanone",
            "cas": "591-78-6",
            "smiles": "CCCCC(=O)C",
            "retention_time": 4.704,
            "response": 4.67,
            "mass_spectrum": None
        },
        {
            "label": "K7_2",
            "group": "ketone",
            "name": "2-heptanone",
            "cas": "110-43-0",
            "smiles": "CCCCCC(=O)C",
            "retention_time": 8.149,
            "response": 3.57,
            "mass_spectrum": None
        },
        {
            "label": "K8_2",
            "group": "ketone",
            "name": "2-octanone",
            "cas": "111-13-7",
            "smiles": "CCCCCCC(=O)C",
            "retention_time": 13.3,
            "response": 2.84,
            "mass_spectrum": None
        },
        {
            "label": "K9_2",
            "group": "ketone",
            "name": "2-nonanone",
            "cas": "821-55-6",
            "smiles": "CCCCCCCC(=O)C",
            "retention_time": 19.517,
            "response": 2.12,
            "mass_spectrum": None
        },

]


def main(root_folder: str, pattern: str):
    specific_folders = glob.glob(root_folder + pattern)
    data = [ca.gc_lc.GCParser.from_Agilent_D_folder(folder) for folder in specific_folders]

    peak_results = []
    for fid, ms in data:
        fid.processor.add(ca.processing.baseline.SectionMinMax())
        ms.processor.add(ca.processing.baseline.SectionMinMax())

    print(data)
    #
    # data = ca.gc_lc.GCMSSignal2D.from_file(folder_ + r"\data.npz")
    # results = ResultGrouper(time_=data.y)
    # for sig in range(len(data)):
    #     results.add_result(sig, process_single(data.get_signal(sig)))
    # TCB = chemistry_lib.find_by_label("TCB")
    # results.set_calibrate(TCB)


def main_first():
    lib_ = ca.gc_lc.GCLibrary('DJW-oxidation')
    lib_.to_JSON(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\new_library.JSON")


if __name__ == "__main__":
    # main_first()
    root_folder = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\standards"
    pattern = r"\DJW-cal-Kn_2-*.D"
    main(root_folder, pattern)
