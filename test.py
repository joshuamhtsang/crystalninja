
import numpy as np
import crystalninja as cn

origin = np.array([0.0,0.0,0.0])


## lattice fcc ${latticeConstant} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
# Ncellx = 3
# Ncelly = 3
# Ncellz = 3
# lattice_constant = 1.0
# orientx = np.array([1.0, 0.0, 0.0])
# orienty = np.array([0.0, 1.0, 0.0])
# orientz = np.array([0.0, 0.0, 1.0])
# pythagorasFactorx = np.linalg.norm(orientx)
# pythagorasFactory = np.linalg.norm(orienty)
# pythagorasFactorz = np.linalg.norm(orientz)
# xdir_celllength = lattice_constant*pythagorasFactorx
# ydir_celllength = lattice_constant*pythagorasFactory
# zdir_celllength = lattice_constant*pythagorasFactorz
# xdir_supercelllength = xdir_celllength*Ncellx
# ydir_supercelllength = ydir_celllength*Ncelly
# zdir_supercelllength = zdir_celllength*Ncellz

## lattice fcc ${latticeConstant} orient x 1 0 0 orient y 0 1 1 orient z 0 -1 1
# Ncellx = 3
# Ncelly = 3
# Ncellz = 3
# lattice_constant = 1.0
# orientx = np.array([ 1.0,  0.0, 0.0])
# orienty = np.array([ 0.0,  1.0, 1.0])
# orientz = np.array([ 0.0, -1.0, 1.0])
# pythagorasFactorx = np.linalg.norm(orientx)
# pythagorasFactory = np.linalg.norm(orienty)
# pythagorasFactorz = np.linalg.norm(orientz)
# xdir_celllength = lattice_constant*pythagorasFactorx
# ydir_celllength = lattice_constant*pythagorasFactory
# zdir_celllength = lattice_constant*pythagorasFactorz
# xdir_supercelllength = xdir_celllength*Ncellx
# ydir_supercelllength = ydir_celllength*Ncelly
# zdir_supercelllength = zdir_celllength*Ncellz

## lattice fcc ${latticeConstant} orient x -1 -1 2 orient y 1 -1 0 orient z 1 1 1
# Ncellx = 3
# Ncelly = 3
# Ncellz = 3
# lattice_constant = 1.0
# orientx = np.array([ -1.0,  -1.0, 2.0])
# orienty = np.array([  1.0,  -1.0, 0.0])
# orientz = np.array([  1.0,  1.0,  1.0])
# pythagorasFactorx = np.linalg.norm(orientx)
# pythagorasFactory = np.linalg.norm(orienty)
# pythagorasFactorz = np.linalg.norm(orientz)
# xdir_celllength = lattice_constant*pythagorasFactorx
# ydir_celllength = lattice_constant*pythagorasFactory
# zdir_celllength = lattice_constant*pythagorasFactorz
# xdir_supercelllength = xdir_celllength*Ncellx
# ydir_supercelllength = ydir_celllength*Ncelly
# zdir_supercelllength = zdir_celllength*Ncellz

## lattice bcc, orientation [1,0,0] [0,1,1] [0,1,1]
Ncellx = 3
Ncelly = 3
Ncellz = 3
lattice_constant = 1.0
orientx = np.array([ 1.0,  1.0, -2.0])
orienty = np.array([ -1.0,  1.0, 0.0])
orientz = np.array([ 1.0,  1.0,  1.0])
pythagorasFactorx = np.linalg.norm(orientx)
pythagorasFactory = np.linalg.norm(orienty)
pythagorasFactorz = np.linalg.norm(orientz)
xdir_celllength = lattice_constant*pythagorasFactorx
ydir_celllength = lattice_constant*pythagorasFactory
zdir_celllength = lattice_constant*pythagorasFactorz
xdir_supercelllength = xdir_celllength*Ncellx
ydir_supercelllength = ydir_celllength*Ncelly
zdir_supercelllength = zdir_celllength*Ncellz




## Use Crystal Ninja.

box1 = cn.box(origin,xdir_supercelllength,ydir_supercelllength,zdir_supercelllength)


#box1.create_crystal("fcc", lattice_constant, orientx, orienty, orientz)
box1.create_crystal("bcc", lattice_constant, orientx, orienty, orientz)

box1.print_lammps_format("out.dat")
