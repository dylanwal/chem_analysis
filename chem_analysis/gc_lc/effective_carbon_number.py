"""

https://doi.org/10.1093/chromsci/23.8.333

Chromatography today
Colin F. Poole and Salwa K. Poole
Pg. 263

"""
import abc

import bigsmiles


class ECNRule(abc.ABC):
    rules = []

    def __init__(self, element: str, value: int | float, description):
        self.element = element
        self.description = description
        self.value = value

    @abc.abstractmethod
    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        ...


class CAliphatic(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'C', 1, 'aliphatic')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        for bond in atom.bonds:
            if bond.bond_order != 1:
                return False
        return True


class COlefinic(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'C', 0.95, '0lefinic')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        for bond in atom.bonds:
            if bond.bond_order == 2 and bond.atom1 == "C" and bond.atom2 == "C":
                return True
        return False


class CAcetylene(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'C', 0.95, 'acetylene')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        for bond in atom.bonds:
            if bond.bond_order == 3 and bond.atom1 == "C" and bond.atom2 == "C":
                return True
        return False


class CCarbonyl(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'C', 0, 'carbonyl')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        for bond in atom.bonds:
            if bond.bond_order == 2 and (bond.atom2 == "O" or bond.atom1 == "O"):
                return True
        return False


class CNitrile(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'C', 0.3, 'nitrile')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        for bond in atom.bonds:
            if bond.bond_order == 3 and (bond.atom2 == "N" or bond.atom1 == "N"):
                return True
        return False


class CAromatic(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'C', 1, 'aromatic')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        if len(atom.root.rings) == 0:
            return False
        # TODO: not correct

        for bond in atom.bonds:
            if bond.bond_order == 2:
                return True

        return False


class OEther(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'O', -1, 'ether')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        if len(atom.bonds) != 2:
            return False
        for bond in atom.bonds:
            if bond.bond_order != 1:
                return False
            if bond.atom2.symbol == "H" or bond.atom1.symbol == "H":
                return False
            if bond.atom1.symbol == "C" and any(bond.bond_order > 1 for bond in bond.atom1.bonds):
                return False
            if bond.atom2.symbol == "C" and any(bond.bond_order > 1 for bond in bond.atom2.bonds):
                return False
        return True


class OPrimaryAlcohol(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'O', -0.5, 'Primary Alcohol')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        if len(atom.bonds) != 2:
            return False

        # check for H and C connection
        H_bond = C_bond = None
        for bond in atom.bonds:
            if bond.atom1 == 'H' or bond.atom2 == 'H':
                H_bond = bond
                continue
            if bond.atom1 == 'C' or bond.atom2 == 'C':
                C_bond = bond

        if H_bond is None or C_bond is None:
            return False

        # check carbon for on two single bonds
        C_atom = C_bond.atom1 if C_bond.atom1 == "C" else C_bond.atom2
        if len(C_atom.bonds) != 2:
            return False
        for bond in C_atom.bonds:
            if bond.bond_order != 1:
                return False

        return True


class OSecondaryAlcohol(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'O', -0.75, 'Secondary Alcohol')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        if len(atom.bonds) != 2:
            return False

        # check for H and C connection
        H_bond = C_bond = None
        for bond in atom.bonds:
            if bond.atom1 == 'H' or bond.atom2 == 'H':
                H_bond = bond
                continue
            if bond.atom1 == 'C' or bond.atom2 == 'C':
                C_bond = bond

        if H_bond is None or C_bond is None:
            return False

        # check carbon for on two single bonds
        C_atom = C_bond.atom1 if C_bond.atom1 == "C" else C_bond.atom2
        if len(C_atom.bonds) != 3:
            return False
        for bond in C_atom.bonds:
            if bond.bond_order != 1:
                return False

        return True


class OTertiaryAlcohol(ECNRule):
    def __init__(self):
        ECNRule.__init__(self, 'O', -0.25, 'Tertiary Alcohol')

    def rule_applies(self, atom: bigsmiles.Atom) -> bool:
        if len(atom.bonds) != 2:
            return False

        # check for H and C connection
        H_bond = C_bond = None
        for bond in atom.bonds:
            if bond.atom1 == 'H' or bond.atom2 == 'H':
                H_bond = bond
                continue
            if bond.atom1 == 'C' or bond.atom2 == 'C':
                C_bond = bond

        if H_bond is None or C_bond is None:
            return False

        # check carbon for on two single bonds
        C_atom = C_bond.atom1 if C_bond.atom1 == "C" else C_bond.atom2
        if len(C_atom.bonds) != 4:
            return False
        for bond in C_atom.bonds:
            if bond.bond_order != 1:
                return False

        return True


RULES = [cls() for cls in ECNRule.__subclasses__()]


def generate_rule_dict():
    dict_ = dict()
    for rule in RULES:
        if rule.element not in dict_:
            dict_[rule.element] = [rule]
        else:
            dict_[rule.element].append(rule)
    return dict_


RULES_dict = generate_rule_dict()


def calculate_effective_carbon_number(smiles: bigsmiles.BigSMILES | str) -> float:
    if isinstance(smiles, str):
        smiles = bigsmiles.BigSMILES(smiles)

    ecn = 0
    for atom in smiles.atoms:
        if atom.symbol in RULES_dict:
            rules = RULES_dict[atom.symbol]  # ensure only correct rules for element is applied
            for rule in rules:
                if rule.rule_applies(atom):
                    ecn += rule.value
                    break
    
    return ecn


def calculate_effective_carbon_number_withlog(smiles: bigsmiles.BigSMILES) -> tuple[float, list[str]]:
    log = []
    ecn = 0
    for atom in smiles.atoms:
        if atom.symbol in RULES_dict:
            rules = RULES_dict[atom.symbol]  # ensure only correct rules for element is applied
            for rule in rules:
                if rule.rule_applies(atom):
                    log.append(type(rule).__name__)
                    ecn += rule.value
                    break

    return ecn, log

"""
N	Primary Amine	-0.50
N	Secondary Amine	-0.75
N	Tertiary Amine	-0.25
Cl	On olefinic C	0.05
Cl	On Alifatic C with 2 Cl atoms	-0.12
"""

if __name__ == '__main__':
    print('calc:', calculate_effective_carbon_number(bigsmiles.BigSMILES("CCCCCCOC(=O)C")), "expected: 7")
    print('calc:', calculate_effective_carbon_number(bigsmiles.BigSMILES("CCCCCCCCC(C)O[Si](C)(C)C")), "expected: 12")
    print('calc:', calculate_effective_carbon_number(bigsmiles.BigSMILES("C[Si](C)(OC(CCCCC(O[Si](C)(C)C)=O)=O)C")), "expected: 10")
    print('calc:', calculate_effective_carbon_number(bigsmiles.BigSMILES("ClC1=C(Cl)C(Cl)=CC=C1")), "expected: 6")
