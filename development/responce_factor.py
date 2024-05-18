
import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\new_library.JSON"
chemistry_lib = ca.mass_spec.GCLibrary.from_JSON(lib_path)


def main(folder_: str):
    data = ca.gc_lc.GCMSSignal2D.from_file(folder_ + r"\data.npz")
    results = ResultGrouper(time_=data.y)
    for sig in range(len(data)):
        results.add_result(sig, process_single(data.get_signal(sig)))
    TCB = chemistry_lib.find_by_label("TCB")
    results.set_calibrate(TCB)


if __name__ == "__main__":
    folder = r"C:\Users\nicep\Desktop\research_wis\data\11\11_23\gc_ms"
    # first_load(folder, "DJW-11-23-*min-TMS.D")
    main(folder)
