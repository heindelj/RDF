import numpy as np
import Frames

class RDF:
    '''
    Class used for calculating the radial distribution function from a collection of frames.
    Also includes functions for plotting the RDF with matplotlib.
    '''
    def __init__(self, frames, nbins=100):
        self.nbins = nbins
        self.frames = frames
        self.frame_by_frame_histogram = np.zeros((len(self.frames, self.nbins)), dtype=np.int32)
        self.histogram = np.zeros(self.nbins, dtype=np.int32)

    def compute_rdf(self, solute_label, solvent_label):
        '''
        Computes total RDF for averaging over all frames.
        '''
        for iFrame in range(len(self.frames)):
            self.compute_single_frame_rdf(iFrame, solute_label, solvent_label)
        

    def compute_single_frame_rdf(self, frame_index, solute_label, solvent_label):
        '''
        Takes two strings, one for the solute element label and the other for the solvent element label.
        These can be the same element for a pair distribution function. Places the histogram into the appropriate
        index of self.frame_by_frame_histogram.
        '''
        try:
            pass
        except IndexError:
            print("Could not access that index. Check how many frames there really are.")