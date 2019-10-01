import pyroshift

chemical_shifts_file = 'PLC_cs.inp'
output_file = 'populations.txt'

#set up the class
chemical_shifts = pyroshift.Isoleucine(chemical_shifts_file,)

#calculate the populations
chemical_shifts.calc_pops_for_all_residues()

#output
chemical_shifts.print_lines()
chemical_shifts.write_lines(output_file)

#write out the latex files!
chemical_shifts.write_latex_pops('latex_populations.tex')
chemical_shifts.write_shift_tex_table('latex_shifts.tex')