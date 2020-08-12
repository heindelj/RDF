#pragma once

#include "Frames.h"
#include <thread>

class RDF
{
public:
	RDF(Frames frames, int nbins, double maxCutoff);

	Frames m_Frames;
	int m_NumBins;
	double m_MaxCutoff;
	double m_BinWidth;

	// parameters
	int numThreads = 16;
	int numSolventAtoms;
	int numFrames;

	std::vector<double> m_Histogram;
	std::vector<std::vector<double>> m_FrameByFrameHistogram;

	void computeRDF(std::string soluteLabel, std::string solventLabel);
	void writeRDF();
private:
	std::vector<std::thread> threads;
	std::vector<double> binRadii;

	void computeSingleFrameRDF(int frameIndex, std::string soluteLabel, std::string solventLabel, double volume);
	double imageDistance(glm::vec3 vec1, glm::vec3 vec2, int frameIndex);
	double smallestAbsoluteComponent(glm::vec3 v);
	std::vector<int> getIndicesOfAtomTypeFromFrame(int frameIndex, std::string& atomLabel);
	void foldImagesIntoBox();
};

