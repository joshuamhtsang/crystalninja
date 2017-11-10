
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

## Problem 5.2.1 in Bulatov and Cai
## lattice bcc, orientation [1,1,-2] [1,1,0] [1,1,1]
#Ncellx = 8
#Ncelly = 19
#Ncellz = 3
#lattice_constant = 2.856
#orientx = np.array([ 1.0,  1.0, -2.0])
#orienty = np.array([ -1.0,  1.0, 0.0])
#orientz = np.array([ 1.0,  1.0,  1.0])
#pythagorasFactorx = np.linalg.norm(orientx)
#pythagorasFactory = np.linalg.norm(orienty)
#pythagorasFactorz = np.linalg.norm(orientz)
#xdir_celllength = lattice_constant*pythagorasFactorx
#ydir_celllength = lattice_constant*pythagorasFactory
#zdir_celllength = lattice_constant*pythagorasFactorz
#xdir_supercelllength = xdir_celllength*Ncellx
#ydir_supercelllength = ydir_celllength*Ncelly
#zdir_supercelllength = zdir_celllength*Ncellz


## Problem 5.2.1 in Bulatov and Cai - adapted for LAMMPS triclinic cell system.
## lattice bcc, orientation [1,1,1] [1,-1,0] [1,1,-2]
Ncellx = 3
Ncelly = 19
Ncellz = 8
lattice_constant = 2.856
orientx = np.array([ 1.0,  1.0, 1.0])
orienty = np.array([ 1.0,  -1.0, 0.0])
orientz = np.array([ 1.0,  1.0,  -2.0])
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
#box1.setTiltFactors(0.0,6.0,0.0)
box1.print_lammps_dump_format("lammpsOut.dat")

## Introduce a screw dislocation dipole.
screw1_pos = np.array([16.3234,43.4188])
screw2_pos = np.array([39.6426,43.4188])
burgersVector = 0.5 * np.array([1,1,1]) * lattice_constant
print "The Burgers Vector is: ", burgersVector
burgersVectorMagnitude = np.linalg.norm(burgersVector)
print "The magnitude of the Burgers Vector is: ", burgersVectorMagnitude
#box1.introduce_screw_dipole_in_xyplane(burgersVectorMagnitude,screw1_pos,screw2_pos)
box1.introduce_screw_dipole_general("y", burgersVectorMagnitude, screw1_pos, screw2_pos)
#box1.setTiltFactors(0.0,0.0,0.0)

box1.print_lammps_dump_format("lammpsOut_screwdipole.dat")

box1.print_lammps_data_format("screwdipole.data")








