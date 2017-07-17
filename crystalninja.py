
import sys
import numpy as np


print "Crystal Ninja - let's slice up some crystal!"





class box:
	def __init__(self, origin_in, x_length_in, y_length_in, z_length_in):
		self.origin   = origin_in
		self.x_length = x_length_in
		self.y_length = y_length_in
		self.z_length = z_length_in
		
		print "The origin is set to: ", self.origin
		print "The x_length is set to: ", self.x_length
		print "The y_length is set to: ", self.y_length
		print "The z_length is set to: ", self.z_length
		
		# Position of the basis atoms in scaled coordinated.
		self.num_basis_atoms = 0
		self.basis = np.array([[0, 0, 0]])
		self.basis = np.delete(self.basis, 0, 0)
		
	
	
	
	def create_crystal(self, structure, a, orientx_in, orienty_in, orientz_in):
		
		print "Creating crystal of structure [", structure, "] within the box."
		
		
		# Set basis atoms for each style.
		if (structure == "fcc"):
			print "Initialising FCC basis atoms..."
			self.num_basis_atoms = 4
			self.basis = np.zeros((self.num_basis_atoms,3))
			self.basis[0] = np.array([0.0, 0.0, 0.0])
			self.basis[1] = np.array([0.5, 0.5, 0.0])
			self.basis[2] = np.array([0.0, 0.5, 0.5])
			self.basis[3] = np.array([0.5, 0.0, 0.5])
			print self.basis
		else:
			print "You haven't specified a valid crystal structure."
			sys.exit(0)
		
		# Normalise the orientation vectors to get the rotation matrix.
		self.orientx  = orientx_in
		self.orienty  = orienty_in
		self.orientz  = orientz_in
		
		print "The orientx is set to: ", self.orientx
		print "The orienty is set to: ", self.orienty
		print "The orientz is set to: ", self.orientz
		
		print "Normalising the orientation vectors and filling rotation matrix..."
		
		self.rotaterow = np.zeros((3,3))
		
		lensq = self.orientx[0]*self.orientx[0] + self.orientx[1]*self.orientx[1] + self.orientx[2]*self.orientx[2]
		length = np.sqrt(lensq)
		self.rotaterow[0,0] = self.orientx[0] / length
		self.rotaterow[0,1] = self.orientx[1] / length
		self.rotaterow[0,2] = self.orientx[2] / length
		
		print self.rotaterow
		
		lensq = self.orienty[0]*self.orienty[0] + self.orienty[1]*self.orienty[1] + self.orienty[2]*self.orienty[2]
		length = np.sqrt(lensq)
		self.rotaterow[1,0] = self.orienty[0] / length
		self.rotaterow[1,1] = self.orienty[1] / length
		self.rotaterow[1,2] = self.orienty[2] / length
		
		print self.rotaterow
		
		lensq = self.orientz[0]*self.orientz[0] + self.orientz[1]*self.orientz[1] + self.orientz[2]*self.orientz[2]
		length = np.sqrt(lensq)
		self.rotaterow[2,0] = self.orientz[0] / length
		self.rotaterow[2,1] = self.orientz[1] / length
		self.rotaterow[2,2] = self.orientz[2] / length
		
		print "The rotation matrix is: \n", self.rotaterow
		
		#self.rotaterow = np.linalg.inv(self.rotaterow)
		
		print "The rotation matrix is: \n", self.rotaterow
		
		#sys.exit(0)
		
		print "Calculating box coordinates of the basis atoms."
		
		self.box_basis = np.zeros((self.num_basis_atoms,3))
		#print self.box_basis
		
		for i in range(0,self.num_basis_atoms):
			self.box_basis[i] = np.matmul(self.rotaterow,self.basis[i])
			#self.box_basis[i] = np.matmul(self.basis[i],self.rotaterow)
		
		print self.box_basis
		
		self.box_basis = np.multiply(self.box_basis, a)
		
		print "The box_basis is:\n", self.box_basis
		
		#sys.exit(0)
		
		print "Filling the box with atoms..."
		
		self.unitcellvector = np.zeros((3,3))
		self.unitcellvector = np.transpose(self.rotaterow) # Holy shit this is important!
		self.unitcellvector = np.multiply(self.unitcellvector,a)
		
		print "The unit cell vector matrix is: \n", self.unitcellvector
		
		#sys.exit(0)
		
		atom_positions = np.array([[0,0,0]])
		atom_positions = np.delete(atom_positions, 0, 0)
		
		box_dimensions_list = [self.x_length, self.y_length, self.z_length]
		longest_box_length = max(box_dimensions_list)
		print "The longest box length is: ", longest_box_length
		max_spawn_index = int( 2*longest_box_length / (a) ) + 1
		print "max_spawn_index = ", max_spawn_index
		
		for i in range(-max_spawn_index,max_spawn_index):
			for j in range(-max_spawn_index,max_spawn_index):
				for k in range(-max_spawn_index,max_spawn_index):
					print "\n(i,j,k) = (", i, ", ", j, ", ", k, ")"
					tmpvec = i*self.unitcellvector[0] + j*self.unitcellvector[1] + k*self.unitcellvector[2]
					print "tmpvec =", tmpvec
					for m in range(0,self.num_basis_atoms):
						tmpvec_atom = np.add(tmpvec,self.box_basis[m])
						new_atom_position = np.array([[ tmpvec_atom[0], tmpvec_atom[1], tmpvec_atom[2] ]]) 
						print "new_atom_position =", new_atom_position
						atom_positions = np.append(atom_positions, new_atom_position, axis=0)
						print "atom_positions: \n", atom_positions
					
					
		
		
		print "atom_positions: \n", atom_positions
		print "There are ", len(atom_positions), " atomic positions generated."
		
		#sys.exit(0)
		
		np.savetxt('data_pretrim.dat', atom_positions)
		
		# Trim away atoms that spawned outside the box dimensions.
		
		deletion_list = []
		
		for i in range(0,len(atom_positions)):
			print "Testing atom with index ", i, " with position: ", atom_positions[i]
			x_position = atom_positions[i,0]
			y_position = atom_positions[i,1]
			z_position = atom_positions[i,2]
			if (x_position > self.x_length or x_position < 0.0):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			elif (y_position > self.y_length or y_position < 0.0):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			elif (z_position > self.z_length or z_position < 0.0):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			else:
				continue
			
		
		print "Deleting ", len(deletion_list), " atoms..."
		print deletion_list
		atom_positions = np.delete(atom_positions,deletion_list,0)
		
		print "There are ", len(atom_positions), " atomic positions left."
		
		np.savetxt('data_posttrim.dat', atom_positions)
		
		print atom_positions
		
		# Trim away atoms that will overlap with existing atoms under PBC.
		
		epsilon = 0.000001
		
		deletion_list = []
		
		for i in range(0,len(atom_positions)):
			print "Testing atom with index ", i, " with position: ", atom_positions[i]
			x_position = atom_positions[i,0]
			y_position = atom_positions[i,1]
			z_position = atom_positions[i,2]
			if (x_position > self.x_length-epsilon):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			elif (y_position > self.y_length-epsilon):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			elif (z_position > self.z_length-epsilon):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			else:
				continue
			
		
		print "Deleting ", len(deletion_list), " atoms..."
		print deletion_list
		atom_positions = np.delete(atom_positions,deletion_list,0)
		
		print "There are ", len(atom_positions), " atomic positions left."
		
		np.savetxt('data_posttrim.dat', atom_positions)
		
		print atom_positions
		
		self.atom_positions = atom_positions
	
	
	
	def print_lammps_format(self, filename):
		print "Printing a LAMMPS-style data file."
		
		ff = open(filename, 'w')
		
		ff.write('ITEM: TIMESTEP\n')
		ff.write('0\n')
		ff.write('ITEM: NUMBER OF ATOMS\n')
		ff.write(str(len(self.atom_positions))+'\n')
		ff.write('ITEM: BOX BOUNDS pp pp pp\n')
		ff.write('0 '+str(self.x_length)+'\n')
		ff.write('0 '+str(self.y_length)+'\n')
		ff.write('0 '+str(self.z_length)+'\n')
		ff.write('ITEM: ATOMS id type x y z\n')
		for i in range(0,len(self.atom_positions)):
			xx = self.atom_positions[i,0]
			yy = self.atom_positions[i,1]
			zz = self.atom_positions[i,2]
			ff.write(str(i)+" "+"1"+" "+str(xx)+" "+str(yy)+" "+str(zz)+"\n")
			
		
		ff.close()
		
		print "Note: you can read in this data file into LAMMPS using: read_data data.txt"
	
	
	
	




