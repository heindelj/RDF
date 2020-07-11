#pragma once

#include "Frames.h"
#include <thread>

class RDF
{
public:
	RDF(Frames frames, int nbins, float maxCutoff);

	Frames m_Frames;
	int m_NumBins;
	double m_MaxCutoff;
	double m_BinWidth = m_MaxCutoff / double(m_NumBins);
	int numThreads = 8;

	std::vector<double> m_Histogram;
	std::vector<std::vector<double>> m_FrameByFrameHistogram;

	void computeRDF(std::string soluteLabel, std::string solventLabel);
	void writeRDF();
private:
	std::vector<std::thread> threads;
	void computeSingleFrameRDF(int frameIndex, std::string soluteLabel, std::string solventLabel, double volume);
	double imageDistance(Vec<double> vec1, Vec<double> vec2, int frameIndex);
	std::vector<int> getIndicesOfAtomTypeFromFrame(int frameIndex, std::string atomLabel);
};

