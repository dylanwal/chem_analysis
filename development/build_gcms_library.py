import json

import chem_analysis as ca


def main():
    chemistry_lib = ca.mass_spec.GCLibrary(name="decane", gc_method="decane")
    file_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\lib.txt"
    with open(file_path, 'r', encoding='UTF-8') as f:
        headers = f.readline()
        headers = headers.split('\t')

        for line in f.readlines():
            line = line.replace('\n', '').split('\t')
            comp = ca.mass_spec.Compound(
                label=line[0],
                group=line[1],
                retention_time=float(line[8]) if line[8] else None,
                smiles=line[4],
                name=line[2],
                cas=line[3],
                response=float(line[9]) if line[9] else None
            )
            chemistry_lib.add_compound(comp)

    lib_dict = chemistry_lib.to_dict(sanitize=True)
    path_ = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.JSON"
    with open(path_, 'w', encoding='UTF-8') as file:
        json.dump(lib_dict, file, indent=4)


def main_load():
    import json
    path_ = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.JSON"
    with open(path_, 'r', encoding='UTF-8') as file:
        lib = ca.mass_spec.GCLibrary.from_dict(json.load(file))

    print(lib)


if __name__ == "__main__":
    # main()
    main_load()
