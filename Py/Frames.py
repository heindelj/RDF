import numpy as np
import sys

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
    
    def __repr__(self):
        return str(self.atomic_symbol + ' ' + np.array2string(self.position).strip('[').strip(']'))

class Frames:
    '''
    Main class for storing frames of molecular data. Contains I/O functions for reading and writing geometries.
    '''
    def __init__(self, xyz_file, a=None, b=None, c=None):
        self.xyz_file = xyz_file
        self.header, self.labels, self.atoms = self.read_geoms(self.xyz_file)
        self.frames = self.get_atoms(self.labels, self.atoms)
        # this is for if the volume is constant as in NVT simulations
        self.cell_parameters = [[a, b, c]]
        self.constant_volume = False
        self.get_cell_parameters()

    def get_cell_parameters(self):
        if None in self.cell_parameters[0]:
            self.cell_parameters = []
            try:
                for header in self.header:
                    self.cell_parameters.append([float(i) for i in header.split()[1:4]])
            except:
                print("Couldn't read box parameters from xyz file. Make sure they are the first things in the comment line, or provide them manually")
                sys.exit(1)
        else:
            self.constant_volume = True
    
    def get_cell_parameters_for_frame(self, frame_index):
        """Gets cell parameters for a particular frame as a list
        """
        if self.constant_volume:
            return self.cell_parameters[0]
        else:
            return self.cell_parameters[frame_index]

    def get_frame_volume(self, frame_index):
        from operator import mul
        from functools import reduce
        return reduce(mul, self.get_cell_parameters_for_frame(frame_index), 1)
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
        return frames
    
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
                line = ifile.readline().split()
                if line and line[0].isdigit():
                    natoms = int(line[0])
                    title = ifile.readline()
                    header.append(str(natoms) + '\n' + ' '.join(line[1:]) + title)
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