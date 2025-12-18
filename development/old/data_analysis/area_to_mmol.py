import csv
import pathlib

import numpy as np


def load_csv(data_path: pathlib.Path):
    data_path = data_path.with_suffix('.csv')
    with open(data_path, mode='r') as file:
        reader = csv.reader(file)
        rows = [row for row in reader]

    # Separate header and data
    times = np.array([float(r) for r in rows[0][slice(1, None)]])
    labels = [r[0] for r in rows[1:]]
    data = [r[slice(1, None)] for r in rows[1:]]
    data = np.array([[float(value) for value in row] for row in data])
    return times, labels, data


def save_csv(filename: pathlib.Path, times: np.ndarray, labels: list[str], mmol: np.ndarray):
    times = np.insert(times, 0, 0)
    text = ",".join([str(i) for i in times]) + "\n"
    for i in range(len(labels)):
        text += ",".join([labels[i]] + [str(i_) for i_ in mmol[i, :]]) + "\n"

    with open(filename, "w") as f:
        f.write(text)


def process_areas(
        data_path: pathlib.Path,
        data_label: str,
        internal_standard: str,
        internal_standard_mmol: float,
        rf: dict[str, float]
):
    times, labels, data = load_csv(data_path / "areas" / data_label)

    # get internal standard
    if internal_standard not in labels:
        raise ValueError(f'Internal standard {internal_standard} not in {labels}')
    internal_standard_index = labels.index(internal_standard)

    rf_array = np.array([rf[label] for label in labels])
    rf_array = rf_array.reshape((len(rf_array), 1))
    mmol = rf_array * (data / data[internal_standard_index, :]) * internal_standard_mmol
    print(mmol)

    # remove IS
    labels.pop(internal_standard_index)
    mmol = np.delete(mmol, internal_standard_index, axis=0)

    save_csv(data_path / "mmols" / (data_label + '.csv'), times, labels, mmol)


def load_ECN(data_path: pathlib.Path):
    data_dict = {}
    with open(data_path, mode='r') as file:
        reader = csv.reader(file)
        for row in reader:
            key = row[0]
            value = int(row[1])
            data_dict[key] = value
    return data_dict


def main_multi():
    data_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\11\11_53")
    internal_standard = "TCB"
    internal_standard_mmol = 0.0058*(15/.1)  # 0.0058 mmol of TCB for .1 ml of rxn ; full rxn is 15 ml

    ECN = load_ECN(pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\ECN_FID_MC.csv"))
    ECN_ref = ECN[internal_standard]
    rf = {label: ECN_ref / ECN_x for label, ECN_x in ECN.items()}

    # get all data_labels and times
    import glob
    import re
    files = glob.glob(str(data_path) + r"\areas\*")
    labels = []
    for file_name in files:
        integers = re.findall(r'\d+', file_name)
        run_num = int(integers[5])
        if run_num in labels:
            raise ValueError(f'duplicate label {run_num}')
        labels.append(run_num)

    # run_xy method on everything
    import multiprocessing

    tasks = [(run_num, data_path, internal_standard, internal_standard_mmol, rf) for run_num in labels]
    with multiprocessing.Pool(multiprocessing.cpu_count()-1) as pool:
        pool.map(process_run, tasks)


def process_run(args):
    run_num, data_path, internal_standard, internal_standard_mmol, rf = args
    data_label = f"DJW-11-53-run_xy{run_num}_fid"
    process_areas(data_path, data_label, internal_standard, internal_standard_mmol, rf)
    print(f"Finished analyzing run_xy {run_num}")


def main():
    data_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\11\11_53")
    data_label = "DJW-11-53-run1_fid"
    internal_standard = "TCB"
    internal_standard_mmol = 0.0058*(15/.1)  # 0.0058 mmol of TCB for .1 ml of rxn ; full rxn is 15 ml

    ECN = load_ECN(pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\ECN_FID_MC.csv"))
    ECN_ref = ECN[internal_standard]
    rf = {label: ECN_ref / ECN_x for label, ECN_x in ECN.items()}

    process_areas(data_path, data_label, internal_standard, internal_standard_mmol, rf)


if __name__ == "__main__":
    # main()
    main_multi()
