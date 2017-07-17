
import numpy as np
import crystalninja as cn

origin = np.array([0.0,0.0,0.0])


## lattice fcc ${latticeConstant} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
#Ncellx = 2
#Ncelly = 2
#Ncellz = 2
#lattice_constant = 1.0
#xdir_celllength = lattice_constant
#ydir_celllength = lattice_constant
#zdir_celllength = lattice_constant
#xdir_supercelllength = xdir_celllength*Ncellx
#ydir_supercelllength = ydir_celllength*Ncelly
#zdir_supercelllength = zdir_celllength*Ncellz
#orientx = np.array([1.0, 0.0, 0.0])
#orienty = np.array([0.0, 1.0, 0.0])
#orientz = np.array([0.0, 0.0, 1.0])

## lattice fcc ${latticeConstant} orient x 1 0 0 orient y 0 1 1 orient z 0 -1 1
#Ncellx = 2
#Ncelly = 2
#Ncellz = 2
#lattice_constant = 1.0
#xdir_celllength = lattice_constant
#ydir_celllength = np.sqrt(2)*lattice_constant
#zdir_celllength = np.sqrt(2)*lattice_constant
#xdir_supercelllength = xdir_celllength*Ncellx
#ydir_supercelllength = ydir_celllength*Ncelly
#zdir_supercelllength = zdir_celllength*Ncellz
#orientx = np.array([ 1.0,  0.0, 0.0])
#orienty = np.array([ 0.0,  1.0, 1.0])
#orientz = np.array([ 0.0, -1.0, 1.0])

## lattice fcc ${latticeConstant} orient x -1 -1 2 orient y 1 -1 0 orient z 1 1 1
Ncellx = 1
Ncelly = 1
Ncellz = 1
lattice_constant = 1.0
xdir_celllength = np.sqrt(6)*lattice_constant
ydir_celllength = np.sqrt(2)*lattice_constant
zdir_celllength = np.sqrt(3)*lattice_constant
xdir_supercelllength = xdir_celllength*Ncellx
ydir_supercelllength = ydir_celllength*Ncelly
zdir_supercelllength = zdir_celllength*Ncellz
orientx = np.array([-1.0, -1.0, 2.0])
orienty = np.array([ 1.0, -1.0, 0.0])
orientz = np.array([ 1.0,  1.0, 1.0])

## Use Crystal Ninja.

box1 = cn.box(origin,xdir_supercelllength,ydir_supercelllength,zdir_supercelllength)
#box1 = cn.box(origin,1,1,1)

box1.create_crystal("fcc", lattice_constant, orientx, orienty, orientz)

box1.print_lammps_format("out.dat")
