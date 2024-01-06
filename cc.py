import os

os.system('ifort -c pnfam_WX.f90')
os.system('make pnfam_nompi.x')
