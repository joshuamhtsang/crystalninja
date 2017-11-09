
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

		# Work out the positional bounds of the 8 corners of the box.
		self.boxbound_xlo = self.origin[0]
		self.boxbound_xhi = self.origin[0] + self.x_length
		self.boxbound_ylo = self.origin[1]
		self.boxbound_yhi = self.origin[1] + self.y_length
		self.boxbound_zlo = self.origin[2]
		self.boxbound_zhi = self.origin[2] + self.z_length

		# Variables for storing the position of the basis atoms in scaled coordinated.
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
		elif (structure == "bcc"):
			self.num_basis_atoms = 2
			self.basis = np.zeros((self.num_basis_atoms,3))
			self.basis[0] = np.array([0.0, 0.0, 0.0])
			self.basis[1] = np.array([0.5, 0.5, 0.5])
		elif (structure == "diamond"):
			print "Initialising Diamond basis atoms..."
			self.num_basis_atoms = 8
			self.basis = np.zeros((self.num_basis_atoms,3))
			self.basis[0] = np.array([0.0, 0.0, 0.0])
			self.basis[1] = np.array([0.0, 0.5, 0.5])
			self.basis[2] = np.array([0.5, 0.0, 0.5])
			self.basis[3] = np.array([0.5, 0.5, 0.0])
			self.basis[4] = np.array([0.25, 0.25, 0.25])
			self.basis[5] = np.array([0.25, 0.75, 0.75])
			self.basis[6] = np.array([0.75, 0.25, 0.75])
			self.basis[7] = np.array([0.75, 0.75, 0.25])
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

		print "Checking that the provides orientx,y,z vectors are sane."
		self.check_orientxyz_orthogonality()
		self.check_orientxyz_righthanded()

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
			print "i = ", i , "/", max_spawn_index
			for j in range(-max_spawn_index,max_spawn_index):
				for k in range(-max_spawn_index,max_spawn_index):
					#print "\n(i,j,k) = (", i, ", ", j, ", ", k, ")"
					tmpvec = i*self.unitcellvector[0] + j*self.unitcellvector[1] + k*self.unitcellvector[2]
					#print "tmpvec =", tmpvec
					# Don't place down atoms if tmpvec is too far from the positive octant.
					if ( tmpvec[0]<(-2.5*a) or tmpvec[1]<(-2.5*a) or tmpvec[2]<(-2.5*a) ):
						continue
					# Don't place down atoms if tmpvec is too far beyond the box boundaries.
					if ( tmpvec[0]>(1.6*self.x_length) or tmpvec[1]>(1.6*self.y_length) or tmpvec[2]>(1.6*self.z_length) ):
						continue
					for m in range(0,self.num_basis_atoms):
						tmpvec_atom = np.add(tmpvec,self.box_basis[m])
						new_atom_position = np.array([[ tmpvec_atom[0], tmpvec_atom[1], tmpvec_atom[2] ]])
						#print "new_atom_position =", new_atom_position
						atom_positions = np.append(atom_positions, new_atom_position, axis=0)
						#print "atom_positions: \n", atom_positions




		print "atom_positions: \n", atom_positions
		print "There are ", len(atom_positions), " atomic positions generated."

		#sys.exit(0)

		np.savetxt('data_pretrim.dat', atom_positions)


		# Trim away atoms that spawned outside the box dimensions.

		atom_positions = self.trim_outsidebox(atom_positions)

		np.savetxt('data_posttrim.dat', atom_positions)


		# Trim away atoms that will overlap with existing atoms under PBC i.e those atoms right on the further out sides of the box.

		atom_positions = self.trim_pbc_overlap(atom_positions)

		np.savetxt('data_posttrim.dat', atom_positions)

		# Shift the position of the pertinent atoms by the user-defined origin.

		atom_positions = self.translate_atom_positions(atom_positions,self.origin)

		self.atom_positions = atom_positions


	# Check orthogonality of the orientx, orienty and orientz vectors.
	def check_orientxyz_orthogonality(self):
		if ( np.dot(self.orientx,self.orienty) != 0 ):
			print "!ERROR! orientx and orienty vectors must be orthogonal, but those provided are not. Please check! :)"
			sys.exit(0)
		if ( np.dot(self.orientx,self.orientz) != 0 ):
			print "!ERROR! orientx and orientz vectors must be orthogonal, but those provided are not. Please check! :)"
			sys.exit(0)
		if ( np.dot(self.orienty,self.orientz) != 0 ):
			print "!ERROR! orienty and orientz vectors must be orthogonal, but those provided are not. Please check! :)"
			sys.exit(0)


	# Check for right-handedness of the 3 orientx,y,z vectors.  x cross y must be in the same direction as z.
	def check_orientxyz_righthanded(self):
		xcrossy = np.cross(self.orientx,self.orienty)
		print xcrossy
		print np.dot(xcrossy,self.orientz)
		if ( np.dot(xcrossy,self.orientz) < 0 ):
			print "!ERROR! The orientx, orienty and orientz vectors do not form a right-handed set.  Please check!  On your right hand, index finger is x, middle finger is y and thumb is z."
			sys.exit(0)



	# Trim away atoms that spawned outside the box dimensions.
	def trim_outsidebox(self, atom_positions):

		deletion_list = []

		epsilon = 0.000001

		for i in range(0,len(atom_positions)):
			#print "Testing atom with index ", i, " with position: ", atom_positions[i]
			x_position = atom_positions[i,0]
			y_position = atom_positions[i,1]
			z_position = atom_positions[i,2]
			if   (x_position > self.x_length+epsilon or x_position < 0.0-epsilon):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			elif (y_position > self.y_length+epsilon or y_position < 0.0-epsilon):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			elif (z_position > self.z_length+epsilon or z_position < 0.0-epsilon):
				#print "Marking atom for deletion..."
				deletion_list = np.append(deletion_list,int(i))
			else:
				continue


		print "Deleting ", len(deletion_list), " atoms..."
		#print deletion_list
		atom_positions = np.delete(atom_positions,deletion_list,0)
		print "There are ", len(atom_positions), " atomic positions left."

		return atom_positions



	# Trim away atoms that will overlap with existing atoms under PBC i.e those atoms right on the further out sides of the box.
	def trim_pbc_overlap(self, atom_positions):

		epsilon = 0.000001

		deletion_list = []

		for i in range(0,len(atom_positions)):
			#print "Testing atom with index ", i, " with position: ", atom_positions[i]
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
		#print deletion_list
		atom_positions = np.delete(atom_positions,deletion_list,0)
		print "There are ", len(atom_positions), " atomic positions left."

		return atom_positions


	def translate_atom_positions(self, atom_positions, tvec):
		print "Translating the atomic positions by [",tvec[0],",",tvec[1],",",tvec[2],"]..."

		atom_positions = atom_positions + tvec

		print "... Translation complete."

		return atom_positions


	def translate_single_atom_position_withpbc(self, x_init, u):

		x_final = np.zeros(3)
		x_final = x_init + u

		if ( x_final[0] > self.boxbound_xhi ):
			x_final[0] -= self.x_length
		if ( x_final[0] < self.boxbound_xlo ):
			x_final[0] += self.x_length
		if ( x_final[1] > self.boxbound_yhi ):
			x_final[1] -= self.y_length
		if ( x_final[1] < self.boxbound_ylo ):
			x_final[1] += self.y_length
		if ( x_final[2] > self.boxbound_zhi ):
			x_final[2] -= self.z_length
		if ( x_final[2] < self.boxbound_zlo ):
			x_final[2] += self.z_length

		return x_final



	def introduce_screw_dipole(self, b, xypos_screw1, xypos_screw2):

		# Check that both dislocations have the same y-coordinate.  This is necessary because the cut between the two dislocations is defined as being parallel to the x-axis.
		if (xypos_screw1[1] != xypos_screw2[1]):
			print "!ERROR! When introducing a screw dislocation dipole, the y-coordinates of both dislocations must be the same."
			sys.exit(0)

		# Work out theta1 and theta2 arrays where:
		#   (1) theta1 : The angles subtended between the positive x-direction and each atom's separation vector R from the dislocation line of screw 1.
		#   (2) theta2 : The same thing as above but for screw 2.
		# Look at Fig. 5.1 in Bulatov and Cai for a diagram.

		num_atoms = len(self.atom_positions)

		theta1 = np.zeros(num_atoms)
		theta2 = np.zeros(num_atoms)

		for i in range(0,num_atoms):
			dx1 = self.atom_positions[i,0] - xypos_screw1[0]
			dy1 = self.atom_positions[i,1] - xypos_screw1[1]
			dx2 = self.atom_positions[i,0] - xypos_screw2[0]
			dy2 = self.atom_positions[i,1] - xypos_screw2[1]
			theta1[i] = self.compute_theta(dx1,dy1)
			theta2[i] = self.compute_theta(dx2,dy2)

		print "theta1 and theta2 arrays successfully computed!"

		# Compute the z-displacements.

		u_z = np.zeros(num_atoms)

		for i in range(0,num_atoms):
			u_z[i] = ( b * (theta1[i]-theta2[i]) ) / (2*np.pi)

		np.savetxt("uz.dat",u_z)
		# Carry out the z-displacements.

		u = np.zeros(3)

		for i in range(0,num_atoms):
			#self.atom_positions[i,2] += u_z[i]
			u = np.array([0,0,u_z[i]])
			self.atom_positions[i] = self.translate_single_atom_position_withpbc(self.atom_positions[i],u)

		print "Screw dipole introduced."




	def compute_theta(self,dx,dy):

		tmp = 0.0

		if ( dy > 0 ):
			tmp = np.arctan2(dy,dx)
		else:
			tmp = np.arctan2(dy,dx) + 2*np.pi

		return tmp



	def print_lammps_dump_format(self, filename):
		print "Printing a LAMMPS-style dump file.  This is useful for viewing by OVITO."

		ff = open(filename, 'w')

		ff.write('ITEM: TIMESTEP\n')
		ff.write('0\n')
		ff.write('ITEM: NUMBER OF ATOMS\n')
		ff.write(str(len(self.atom_positions))+'\n')
		ff.write('ITEM: BOX BOUNDS pp pp pp\n')
		ff.write(str(self.boxbound_xlo)+' '+str(self.boxbound_xhi)+'\n')
		ff.write(str(self.boxbound_ylo)+' '+str(self.boxbound_yhi)+'\n')
		ff.write(str(self.boxbound_zlo)+' '+str(self.boxbound_zhi)+'\n')
		ff.write('ITEM: ATOMS id type x y z\n')
		for i in range(0,len(self.atom_positions)):
			xx = self.atom_positions[i,0]
			yy = self.atom_positions[i,1]
			zz = self.atom_positions[i,2]
			ff.write(str(i+1)+" "+"1"+" "+str(xx)+" "+str(yy)+" "+str(zz)+"\n")


		ff.close()


	def print_lammps_data_format(self, filename):
		print "Printing a LAMMPS-style data file.  This file can be read by LAMMPS for an initial configuration using the 'read_data' LAMMPS command."

		ff = open(filename, 'w')

		ff.write('MESSAGE FROM CRYSTAL NINJA: INSERT YOUR OWN COMMENT LINE HERE.\n')
		ff.write('\n')
		ff.write(str(len(self.atom_positions))+' atoms\n')
		ff.write('1 atom types\n')
		ff.write(str(self.boxbound_xlo)+' '+str(self.boxbound_xhi)+'xlo xhi'+'\n')
		ff.write(str(self.boxbound_ylo)+' '+str(self.boxbound_yhi)+'ylo yhi'+'\n')
		ff.write(str(self.boxbound_zlo)+' '+str(self.boxbound_zhi)+'zlo zhi'+'\n')
		ff.write('\n')
		ff.write(' Atoms\n')
		ff.write('\n')
		for i in range(0,len(self.atom_positions)):
			xx = self.atom_positions[i,0]
			yy = self.atom_positions[i,1]
			zz = self.atom_positions[i,2]
			ff.write(str(i+1)+" "+"1"+" "+str(xx)+" "+str(yy)+" "+str(zz)+"\n")


		ff.close()
