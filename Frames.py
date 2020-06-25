import numpy as np

class Atom:
    '''
    Atom class. Stores atom identity and position.
    '''
    def __init__(self, atomic_symbol, position):
        '''
        atomic_symbol is a string for the element symbol.
        position is a 3D numpy array of cartesian coordinates
        '''
        self.atomic_symbol = atomic_symbol
        self.position = np.array(position)

class Frames:
    '''
    Main class for storing frames of molecular data. Contains I/O functions for reading and writing geometries.
    '''
    def __init__(self, xyz_file):
        self._xyz_file = xyz_file
        self.header, self.labels, self.atoms = self.read_geoms(self._xyz_file)
        self.frames = self.get_atoms(self.labels, self.atoms)
    
    def __repr__(self):
        return self.frames
    
    def get_number_of_atom_type_from_frame(self, frame_index, atom_label):
        '''
        Get total number of atoms of a particular type.
        '''
        num_atoms = 0
        for atom in self.frames[frame_index]:
            if atom_label == atom.atomic_symbol:
                num_atoms += 1
        if num_atoms > 0:
            return num_atoms
        return None
    
    def get_total_number_of_atoms_from_frame(self, frame_index):
        '''
        Get total number of atoms from a frame.
        '''
        return len(self.frames[frame_index])

    
    @staticmethod
    def get_atoms(labels, atoms):
        '''
        Takes a list of atom labels and the coordinates of those atoms.
        labels: list of lists of strings
        atoms: list of lists of numpy arrays
        frames: list of lists of Atom objects
        '''
        frames = []
        for iFrame, frame_labels in enumerate(labels):
            frame = []
            for iAtom, atom_label in enumerate(frame_labels):
                frame.append(Atom(atom_label, atoms[iFrame][iAtom]))
            frames.append(frame)
    
    @staticmethod
    def read_geoms(geom):
        '''
        Reads a file containing a large number of XYZ formatted files concatenated together and splits them
        into an array of arrays of vectors (MxNx3) where M is the number of geometries, N is the number of atoms,
        and 3 is from each x, y, z coordinate.
        '''
        allCoords = []
        atomLabels = []
        header = []
        with open(geom) as ifile:
            while True:
                atomLabels__ = []
                line = ifile.readline()
                if line.strip().isdigit():
                    natoms = int(line)
                    title = ifile.readline()
                    header.append(str(natoms) + '\n' + title)
                    coords = np.zeros([natoms, 3], dtype="float64")
                    for x in coords:
                        line = ifile.readline().split()
                        atomLabels__.append(line[0])
                        x[:] = list(map(float, line[1:4]))
                    allCoords.append(coords)
                    atomLabels.append(atomLabels__)
                if not line:
                    break
        return header, atomLabels, allCoords

    def write_geoms(self, ofile=None):
        """
        args: header, labels, and coords as output by read_geoms
        return: no return

        Writes the the moelcules to stdout in xyz format if no ofile is specified.
        Otherwise, write the geometries to ofile.
        """
        for i, head in enumerate(self.header):
            output = head
            for j, vec in enumerate(self.atoms[i]):
                output += (self.labels[i][j] + " " +
                        np.array2string(np.array(vec), precision=14, separator=' ', suppress_small=True).strip('[]') + '\n')
            if ofile is None:
                print(output)
            elif i == 0:
                with open(ofile, 'w') as f:
                    f.write(output)
            else:
                with open(ofile, 'a') as f:
                    f.write(output)