import os
from ase.io import gen, vasp

##############
# structures #
##############
path = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))
mef = vasp.read_vasp(path + "/reference_files/CONTCAR_mef")
cf4 = vasp.read_vasp(path + "/reference_files/CONTCAR_cf4")
amorphous = gen.read_gen(path + "/reference_files/amorphous_base.gen")
xtl_n = vasp.read_vasp(path + "/reference_files/CONTCAR_nrich")
xtl_si = vasp.read_vasp(path + "/reference_files/CONTCAR_sirich")
xtl2x2 = gen.read_gen(path + "/reference_files/2x2xtl.gen")
xtl2x2_sifterm = gen.read_gen(path + "/reference_files/2x2xtl_sifterm.gen")
heavy_bomb = vasp.read_vasp(path + "/reference_files/CONTCAR_heavy_bombard")
bulk222 = vasp.read_vasp(path + "/reference_files/CONTCAR_222bulk")
annealed = vasp.read_vasp(path + "/reference_files/CONTCAR_annealed_unitcell")

