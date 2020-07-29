from Frames import Frames
from RDF import RDF

if __name__ == '__main__':
    trajectory = Frames("C:/Users/jhein/OneDrive/Documents/Research/Sotiris/BSSE_MBE/RDF_Code/test_trajectory/traj_short.xyz")
    OO_rdf = RDF(trajectory, nbins=150, max_cutoff=10.0)
    OO_rdf.compute_rdf('O', 'O')
    OO_rdf.plot()