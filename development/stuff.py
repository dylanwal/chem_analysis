import rxn
print(rxn.__file__)

print(help(rxn))
table = rxn.data_structure
print(help(table))
print(table.PERIODIC_TABLE.print_table())
element = table.PERIODIC_TABLE[10]
print(element.atomic_mass)
# print(element.isotope_most_abundant_mass())
print(help(element))

result = rxn.calculators.find_combinations([("C",0, 10),("O",0,10),("H",0, 10)], 100, 0.5)
print(result)