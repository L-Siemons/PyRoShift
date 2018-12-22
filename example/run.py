import pyroshift

input = pyroshift.Input('PLC_cs.inp')
chemical_shifts = pyroshift.isoleucine('PLC_cs.inp')
print chemical_shifts.calc_pops_for_all_residues()
