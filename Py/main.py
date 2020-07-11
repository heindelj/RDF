from Frames import Frames
from RDF import RDF

if __name__ == '__main__':
    trajectory = Frames("test_trajectory/traj_short.xyz")
    OO_rdf = RDF(trajectory, nbins=150)
    OO_rdf.compute_rdf('O', 'O')
    print(OO_rdf.histogram)
    OO_rdf.plot()