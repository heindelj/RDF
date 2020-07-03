from Frames import Frames
from RDF import RDF

if __name__ == '__main__':
    trajectory = Frames("test_trajectory/traj.xyz")
    OO_rdf = RDF(trajectory, nbins=150)
    OO_rdf.compute_rdf('O', 'O')
    OO_rdf.plot()