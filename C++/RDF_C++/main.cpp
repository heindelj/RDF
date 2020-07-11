#include "Frames.h"
#include "RDF.h"
#include "Timer.h"

int main()
{
	Timer timer;
	timer.begin();
	Frames trajectory = Frames("C:\\Users\\jhein\\OneDrive\\Documents\\Research\\Sotiris\\BSSE_MBE\\RDF_Code\\test_trajectory\\traj.xyz");
	std::cout << timer.end() / 1000000 << " seconds to read input file." << std::endl;
	RDF OO_rdf = RDF(trajectory, 150, 12.0);

	OO_rdf.numThreads = 32; // gets about 110 frames per second

	timer.begin();
	OO_rdf.computeRDF("O", "O");
	std::cout << float(timer.end()) / 1000000.0 << " seconds to compute RDF." << std::endl;
	std::cout << trajectory.m_Frames.size() / (timer.end() / 1000000.0) << " frames per second." << std::endl;
	OO_rdf.writeRDF();
	std::cin.get();
}