from itertools import product

# Define a dictionary of elements with their atomic weights (in g/mol)
element_weights = {
    'H': 1.008,
    # 'He': 4.002602,
    # 'Li': 6.94,
    # 'Be': 9.0122,
    # 'B': 10.81,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    # 'F': 18.998,
    # 'Ne': 20.180,
    # 'Na': 22.990,
    # 'Mg': 24.305,
    # 'Al': 26.982,
    # 'Si': 28.085,
    # 'P': 30.974,
    # 'S': 32.06,
    # 'Cl': 35.45,
    # 'Ar': 39.948,
    # 'K': 39.098,
    # 'Ca': 40.078,
    # 'Sc': 44.956,
    # 'Ti': 47.867,
    # 'V': 50.9415,
    # 'Cr': 51.9961,
    # 'Mn': 54.938044,
    # 'Fe': 55.845,
    # 'Co': 58.933194,
    # 'Ni': 58.933,
    # 'Cu': 63.546,
    # 'Zn': 65.38,
    # 'Ga': 69.723,
    # 'Ge': 72.630,
    # 'As': 74.921595,
    # 'Se': 78.971,
    # 'Br': 79.904,
    # 'Kr': 83.798,
    # 'Rb': 85.4678,
    # 'Sr': 87.621,
    # 'Y': 88.90584,
    # 'Zr': 91.224,
    # 'Nb': 92.90637,
    # 'Mo': 95.951,
    # 'Tc': 98,
    # 'Ru': 101.07,
    # 'Rh': 102.90550,
    # 'Pd': 106.42,
    # 'Ag': 107.8682,
    # 'Cd': 112.414,
    # 'In': 114.818,
    # 'Sn': 118.710,
    # 'Sb': 121.760,
    # 'Te': 127.60,
    # 'I': 126.90447,
    # 'Xe': 131.293,
    # 'Cs': 132.90545196,
    # 'Ba': 137.327,
    # 'La': 138.90547,
    # 'Ce': 140.116,
    # 'Pr': 140.90766,
    # 'Nd': 144.242,
    # 'Pm': 145,
    # 'Sm': 150.36,
    # 'Eu': 151.964,
    # 'Gd': 157.25,
    # 'Tb': 158.92535,
    # 'Dy': 162.500,
    # 'Ho': 164.93033,
    # 'Er': 167.259,
    # 'Tm': 168.93422,
    # 'Yb': 173.04,
    # 'Lu': 174.9668,
    # 'Hf': 178.49,
    # 'Ta': 180.94788,
    # 'W': 183.84,
    # 'Re': 186.207,
    # 'Os': 190.23,
    # 'Ir': 192.217,
    # 'Pt': 195.084,
    # 'Au': 196.966569,
    # 'Hg': 200.592,
    # 'Tl': 204.38,
    # 'Pb': 207.2,
    # 'Bi': 208.98040,
    # 'Po': 209,
    # 'At': 210,
    # 'Rn': 222,
    # 'Fr': 223,
    # 'Ra': 226,
    # 'Ac': 227,
    # 'Th': 232.0377,
    # 'Pa': 231.03588,
    # 'U': 238.02891,
    # 'Np': 237,
    # 'Pu': 244,
    # 'Am': 243,
    # 'Cm': 247,
    # 'Bk': 247,
    # 'Cf': 251,
    # 'Es': 252,
    # 'Fm': 257,
    # 'Md': 258,
    # 'No': 259,
    # 'Lr': 262,
    # 'Rf': 267,
    # 'Db': 270,
    # 'Sg': 271,
    # 'Bh': 270,
    # 'Hs': 277,
    # 'Mt': 276,
    # 'Ds': 281,
    # 'Rg': 282,
    # 'Cn': 285,
    # 'Nh': 286,
    # 'Fl': 289,
    # 'Mc': 288,
    # 'Lv': 293,
    # 'Ts': 294,
    # 'Og': 294,
}

def find_combinations(target_weight, tolerance, max_elements=15):
    combinations = []

    # Generate possible counts for each element (1 to max_elements)
    for counts in product(range(max_elements + 1), repeat=len(element_weights)):
        total_weight = sum(count * weight for count, weight in zip(counts, element_weights.values()))

        # Check if the total weight is within the target weight +/- tolerance
        if target_weight - tolerance <= total_weight <= target_weight + tolerance:
            # Create a readable format of the combination
            combination = {element: count for element, count in zip(element_weights.keys(), counts) if count > 0}
            combinations.append((combination, total_weight))

    return combinations


from collections import defaultdict


def find_combinations_greedy(target_weight, tolerance, max_elements=15):
    combinations = []

    # Sort elements by atomic weight
    sorted_elements = sorted(element_weights.items(), key=lambda x: x[1])

    # Greedily try to find combinations
    for count in range(1, max_elements + 1):
        for element_combination in product(sorted_elements, repeat=count):
            element_dict = defaultdict(int)
            total_weight = 0

            for element, _ in element_combination:
                element_dict[element] += 1
                total_weight += element_weights[element]

            # Check if the total weight is within the target weight +/- tolerance
            if target_weight - tolerance <= total_weight <= target_weight + tolerance:
                combinations.append((dict(element_dict), total_weight))

    return combinations

def find_combinations_trimming(target_weight, tolerance, max_elements=15):
    combinations = []

    # Generate possible counts for each element (1 to max_elements)
    for counts in product(range(max_elements + 1), repeat=len(element_weights)):
        total_weight = sum(count * weight for count, weight in zip(counts, element_weights.values()))

        # Check if the total weight exceeds the target weight + tolerance
        if total_weight > target_weight + tolerance:
            continue  # Skip this combination

        # Check if the total weight is within the target weight +/- tolerance
        if target_weight - tolerance <= total_weight <= target_weight + tolerance:
            # Create a readable format of the combination
            combination = {element: count for element, count in zip(element_weights.keys(), counts) if count > 0}
            combinations.append((combination, total_weight))

    return combinations

max_bonds = {
    'C': 4,  # Carbon
    'H': 1,  # Hydrogen
    'N': 3,  # Nitrogen
    'O': 2   # Oxygen
}

def find_combinations_2(target_weight, tolerance):
    combinations = []

    # Loop through possible counts for each element
    for c_count in range(max_bonds['C'] + 1):
        for h_count in range(max_bonds['H'] + 1):
            for n_count in range(max_bonds['N'] + 1):
                for o_count in range(max_bonds['O'] + 1):
                    # Calculate total bonds
                    total_bonds = (c_count * max_bonds['C'] +
                                   h_count * max_bonds['H'] +
                                   n_count * max_bonds['N'] +
                                   o_count * max_bonds['O'])

                    # Check if total bonds exceed available bonds
                    if total_bonds > (c_count * max_bonds['C'] + h_count * max_bonds['H'] +
                                      n_count * max_bonds['N'] + o_count * max_bonds['O']):
                        continue  # Skip this combination

                    # Calculate total weight
                    total_weight = (c_count * element_weights['C'] +
                                    h_count * element_weights['H'] +
                                    n_count * element_weights['N'] +
                                    o_count * element_weights['O'])

                    # Check if the total weight is within the target weight +/- tolerance
                    if target_weight - tolerance <= total_weight <= target_weight + tolerance:
                        combination = {
                            'C': c_count,
                            'H': h_count,
                            'N': n_count,
                            'O': o_count
                        }
                        combinations.append((combination, total_weight))

    return combinations

# Example usage
target_molecular_weight = 67.0  # target molecular weight in g/mol
tolerance = 2.0  # tolerance in g/mol

result = find_combinations(target_molecular_weight, tolerance)
print(len(result))
# for combo, weight in result:
#     print(f"Combination: {combo}, Total Weight: {weight:.3f} g/mol")

result = find_combinations_trimming(target_molecular_weight, tolerance)
print(len(result))
# for combo, weight in result:
#     print(f"Combination: {combo}, Total Weight: {weight:.3f} g/mol")


result = find_combinations_2(target_molecular_weight, tolerance)
print(len(result))
for combo, weight in result:
    print(f"Combination: {combo}, Total Weight: {weight:.3f} g/mol")