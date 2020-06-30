import numpy as np

class RDF:
    '''
    Class used for calculating the radial distribution function from a list of frames (i.e. Frame.frames).
    Must provide the lattice parameters of the cell.
    Also includes functions for plotting the RDF with matplotlib.
    '''
    def __init__(self, frames, a=1.0, b=1.0, c=1.0, nbins=100, max_cutoff=12.0):
        self.nbins = nbins
        self.frames = frames
        self.a = a
        self.b = b
        self.c = c
        self.volume = self.a * self.b * self.c
        self.max_cutoff = max_cutoff
        self.bin_width = self.max_cutoff / self.nbins
        self.frame_by_frame_histogram = np.zeros((len(self.frames, self.nbins)), dtype=np.int32)
        self.histogram = np.zeros(self.nbins, dtype=np.int32)

    def compute_rdf(self, solute_label, solvent_label):
        '''
        Computes total RDF for averaging over all frames.
        '''
        for iFrame in range(len(self.frames)):
            self.compute_single_frame_rdf(iFrame, solute_label, solvent_label)
        #some kind of np.mean thing here
        

    def compute_single_frame_rdf(self, frame_index, solute_label, solvent_label):
        '''
        Takes two strings, one for the solute element label and the other for the solvent element label.
        These can be the same element for a pair distribution function. Places the histogram into the appropriate
        index of self.frame_by_frame_histogram.
        '''
        num_atoms_of_type = len(self.get_indices_of_atom_type_from_frame(frame_index, solvent_label))
        try:
            for iSolute in self.get_indices_of_atom_type_from_frame(frame_index, solute_label):
                for iSolvent in self.get_indices_of_atom_type_from_frame(frame_index, solvent_label):
                    if iSolute != iSolvent: # need to account for PBC here in some way
                        r = np.linalg.norm(self.frames[frame_index][iSolute].position - self.frames[frame_index][iSolvent].position)
                        bin_index = int(r / self.bin_width)
                        self.frame_by_frame_histogram[frame_index][bin_index] += 1.0
            for bin_index in range(self.nbins):
                self.frame_by_frame_histogram[frame_index][bin_index] /= \
                    ((num_atoms_of_type / self.volume * 4 * np.pi * ((bin_index + 1) * self.bin_width )**2) )

        except IndexError:
            print("Could not access that index. Check how many frames there really are.")
    
    def get_indices_of_atom_type_from_frame(self, frame_index, atom_label):
        '''
        Get indices of atoms of a particular type and return as list.
        '''
        atom_indices = []
        for iAtom, atom in enumerate(self.frames[frame_index]):
            if atom_label == atom.atomic_symbol:
                atom_indices.append(iAtom)
        if len(atom_indices) > 0:
            return atom_indices
        return None
    
    def get_total_number_of_atoms_from_frame(self, frame_index):
        '''
        Get total number of atoms from a frame.
        '''
        return len(self.frames[frame_index])