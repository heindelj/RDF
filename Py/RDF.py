import numpy as np
import multiprocessing as mp
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import tqdm

class RDF:
    '''
    Class used for calculating the radial distribution function from a Frames object.
    Also includes functions for plotting the RDF with matplotlib.
    '''
    def __init__(self, frames, nbins=100, max_cutoff=12.0):
        self.nbins = nbins
        self.Frames = frames
        self.max_cutoff = max_cutoff
        self.bin_width = self.max_cutoff / self.nbins
        self.frame_by_frame_histogram = np.zeros((len(self.Frames.frames), self.nbins), dtype=np.float32)
        self.histogram = None

    def compute_rdf(self, solute_label, solvent_label):
        """Computes total rdf averaging over all bins and stores result in self.histogram

        Args:
            solute_label (String): Element symbol for solute (X)
            solvent_label (String): Element symbol for solvent
        """
        #num_cores = multiprocessing.cpu_count()
        num_frames = len(self.Frames.frames)
        print("Found " + str(num_frames) + " frames:")

        # get arguments for process parallel calculation of RDF
        #args__ = []
        #for iFrame in range(num_frames):
        #    args__.append([iFrame, solute_label, solvent_label, self.Frames.get_frame_volume(iFrame)])
        #args = tqdm.tqdm(args__)
        #with mp.Pool(processes=mp.cpu_count()) as pool:
        #    pool.starmap(self.compute_single_frame_rdf, args)

        for iFrame in tqdm.tqdm(range(num_frames)):
            self.compute_single_frame_rdf(iFrame, solute_label, solvent_label, self.Frames.get_frame_volume(iFrame))
            if iFrame > 0 and (iFrame % 5000 == 0):
                self.histogram = np.mean(self.frame_by_frame_histogram[:iFrame], axis=0)
                self.plot()
        self.histogram = np.mean(self.frame_by_frame_histogram, axis=0)
        

    def compute_single_frame_rdf(self, frame_index, solute_label, solvent_label, volume):
        '''
        Takes two strings, one for the solute element label and the other for the solvent element label.
        These can be the same element for a pair distribution function. Places the histogram into the appropriate
        index of self.frame_by_frame_histogram.
        '''
        num_atoms_of_type = len(self.get_indices_of_atom_type_from_frame(frame_index, solvent_label))
        try:
            # do this branching to avoid double-checking pairs when solute and solvent are the same
            solute_indices = self.get_indices_of_atom_type_from_frame(frame_index, solute_label)
            if solute_label == solvent_label:
                solvent_indices = solute_indices
                for i, iSolute in enumerate(solute_indices[0:-1]):
                    for iSolvent in solvent_indices[i+1:]:
                        r = self.image_distance(self.Frames.frames[frame_index][iSolute].position \
                                                ,self.Frames.frames[frame_index][iSolvent].position, frame_index)
                        if r < self.max_cutoff:
                            bin_index = int(r / self.bin_width)
                            self.frame_by_frame_histogram[frame_index][bin_index] += 2.0
            else:
                solvent_indices = self.get_indices_of_atom_type_from_frame(frame_index, solvent_label)
                for iSolute in solute_indices:
                    for iSolvent in solvent_indices:
                        r = self.image_distance(self.Frames.frames[frame_index][iSolute].position \
                                                ,self.Frames.frames[frame_index][iSolvent].position, frame_index)
                        if r < self.max_cutoff:
                            bin_index = int(r / self.bin_width)
                            self.frame_by_frame_histogram[frame_index][bin_index] += 1.0
            # renormalize bins to ideal gas probability distribution
            for bin_index in range(self.nbins):
                self.frame_by_frame_histogram[frame_index][bin_index] /= \
                    (4 * np.pi * (num_atoms_of_type * (num_atoms_of_type - 1) / volume) * ((bin_index + 1) * self.bin_width )**2 * self.bin_width)

        except IndexError:
            print("Could not access that index. Check how many frames there really are.")
    
    def image_distance(self, vec1, vec2, frame_index):
        """Returns distance to nearest image of a molecule when PBC are used

        Args:
            vec1 (3d numpy array or list): Position of one atom
            vec2 (3d numpy array or list): Position of second atom
        """
        def closest_image(side_length, distance_component):
            """Finds the distance in x, y, or z to the nearest image.

            Args:
                side_length (float): box side length
                distance_component (float): x, y, or
                 z component of distance vector

            Returns:
                float: component of distance vector to nearest image
            """
            test_distance = side_length * 0.5
            if abs(distance_component) <= test_distance:
                return distance_component
            if distance_component < (-test_distance):
                return distance_component + int(abs((distance_component - test_distance)/ side_length)) * side_length
            if distance_component > (test_distance):
                return distance_component - int(abs((distance_component + test_distance)/ side_length)) * side_length
            

        #should be vectorized but doesn't make any performance differece
        distance_vector = vec1 - vec2
        cell_parameters = self.Frames.get_cell_parameters_for_frame(frame_index)
        dx = closest_image(cell_parameters[0], distance_vector[0])
        dy = closest_image(cell_parameters[1], distance_vector[1])
        dz = closest_image(cell_parameters[2], distance_vector[2])
        return np.sqrt(dx**2 + dy**2 + dz**2)
        #return np.linalg.norm(vec1 - vec2)

    def plot(self):
        """
        Plots the RDF using matplotlib.
        """
        ax1 = plt.figure().gca()
        ax1.plot(np.linspace(0.0, self.max_cutoff, num=self.nbins), self.histogram)
        plt.xlabel("Pair Distance ($\AA$)")
        plt.ylabel("g(r)")
        plt.show()

    def get_indices_of_atom_type_from_frame(self, frame_index, atom_label):
        """gets indices of atom type from a particular frame

        Args:
            frame_index (int): frame index
            atom_label (string): element symbol being looked for

        Returns:
            [list of ints]: all indices to atom of requested type
        """
        atom_indices = []
        for iAtom, atom in enumerate(self.Frames.frames[frame_index]):
            if atom_label == atom.atomic_symbol:
                atom_indices.append(iAtom)
        if len(atom_indices) > 0:
            return atom_indices
        return None
    
    def get_total_number_of_atoms_from_frame(self, frame_index):
        '''
        Get total number of atoms from a frame.
        '''
        return len(self.Frames.frames[frame_index])