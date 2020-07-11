#include "RDF.h"
#include <cstdlib>     // abs
#define _USE_MATH_DEFINES
#include <math.h>      // sqrt, pi

RDF::RDF(Frames frames, int nbins, float maxCutoff) :
	m_Frames(frames), m_NumBins(nbins), m_MaxCutoff(maxCutoff)
{
	std::vector<double> histogram(m_NumBins);
	m_Histogram = histogram;

	// fill in the frame_by_frame_histrogram
	for (int i = 0; i < m_Frames.m_Frames.size(); i++)
		m_FrameByFrameHistogram.push_back(m_Histogram);

	m_BinWidth = m_MaxCutoff / float(m_NumBins);
}

void RDF::writeRDF()
{
	for (int i = 0; i < m_Histogram.size(); i++)
		std::cout << m_BinWidth * (double(i)+1.0) << "   " << m_Histogram[i] << std::endl;
}

void RDF::computeRDF(std::string soluteLabel, std::string solventLabel)
{
	int numThreadedBlocks = int(m_Frames.m_Frames.size() / numThreads);
	for (int iBlock = 0; iBlock < numThreadedBlocks; iBlock++)
	{
		for (int iThread = 0; iThread < numThreads; iThread++)
			threads.emplace_back([=]() {computeSingleFrameRDF(iBlock * numThreads + iThread, soluteLabel, solventLabel, m_Frames.getFrameVolume(iBlock * numThreads + iThread)); });

		for (auto& thread : threads)
			thread.join();

		threads.clear();
	}

	// finish up the ones which didn't fit in a block
	for (int frameIndex = numThreadedBlocks * numThreads; frameIndex < m_Frames.m_Frames.size(); frameIndex++)
		computeSingleFrameRDF(frameIndex, soluteLabel, solventLabel, m_Frames.getFrameVolume(frameIndex));

	for (int iFrame = 0; iFrame < m_FrameByFrameHistogram.size(); iFrame++)
	{
		for (int iBin = 0; iBin < m_Histogram.size(); iBin++)
			m_Histogram[iBin] += m_FrameByFrameHistogram[iFrame][iBin] / m_FrameByFrameHistogram.size();
	}
}

void RDF::computeSingleFrameRDF(int frameIndex, std::string soluteLabel, std::string solventLabel, double volume)
{
	if ((frameIndex % 10000) == 0)
		std::cout << "Computing Frame Number: " << frameIndex << std::endl;
	std::vector<int> soluteIndices = getIndicesOfAtomTypeFromFrame(frameIndex, soluteLabel);
	std::vector<int> solventIndices = getIndicesOfAtomTypeFromFrame(frameIndex, solventLabel);

	if (soluteLabel == solventLabel)
	{
		for (int iSolute = 0; iSolute < soluteIndices.size() - 1; iSolute++)
		{
			// figure out if we're going to use the pair list or not
			for (int iSolvent = iSolute + 1; iSolvent < solventIndices.size(); iSolvent++)
			{
				double r = imageDistance(m_Frames.m_Frames[frameIndex][soluteIndices[iSolute]].m_Position, m_Frames.m_Frames[frameIndex][solventIndices[iSolvent]].m_Position, frameIndex);
				if (r < m_MaxCutoff && r > 0.0)
				{
					int binIndex = int(r / m_BinWidth);
					m_FrameByFrameHistogram[frameIndex][binIndex] += 1.0;
				}
			}
		}
	}
	else
	{
		for (int iSolute : soluteIndices)
		{
			for (int iSolvent : solventIndices)
			{
				double r = imageDistance(m_Frames.m_Frames[frameIndex][iSolute].m_Position, m_Frames.m_Frames[frameIndex][iSolvent].m_Position, frameIndex);
				if (r < m_MaxCutoff && r > 0.0)
				{
					int binIndex = int(r / m_BinWidth);
					m_FrameByFrameHistogram[frameIndex][binIndex] += 1.0;
				}
			}
		}
	}

	// renormalize to ideal gas RDF
	int numAtomsOfType = getIndicesOfAtomTypeFromFrame(frameIndex, solventLabel).size();
	for (int binIndex = 0; binIndex < m_NumBins; binIndex++)
		m_FrameByFrameHistogram[frameIndex][binIndex] /= (4 * M_PI * double(numAtomsOfType) / volume
			* double((binIndex + 1) * (binIndex + 1)) * (m_BinWidth * m_BinWidth)) * m_BinWidth;
}

double RDF::imageDistance(Vec<double> vec1, Vec<double> vec2, int frameIndex)
{
	Vec<double> distanceVector = vec1 - vec2;
	// branchless calculation of components
	double dx = (distanceVector[0]) * (abs(distanceVector[0]) <= (m_Frames.m_CellParameters[frameIndex][0] / 2.0))
		+ (m_Frames.m_CellParameters[frameIndex][0] - abs(distanceVector[0])) * (abs(distanceVector[0]) > (m_Frames.m_CellParameters[frameIndex][0] / 2.0));
	double dy = (distanceVector[1]) * (abs(distanceVector[1]) <= (m_Frames.m_CellParameters[frameIndex][1] / 2.0))
		+ (m_Frames.m_CellParameters[frameIndex][1] - abs(distanceVector[1])) * (abs(distanceVector[1]) > (m_Frames.m_CellParameters[frameIndex][1] / 2.0));
	double dz = (distanceVector[2]) * (abs(distanceVector[2]) <= (m_Frames.m_CellParameters[frameIndex][2] / 2.0))
		+ (m_Frames.m_CellParameters[frameIndex][2] - abs(distanceVector[2])) * (abs(distanceVector[2]) > (m_Frames.m_CellParameters[frameIndex][2] / 2.0));
	return sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<int> RDF::getIndicesOfAtomTypeFromFrame(int frameIndex, std::string atomLabel)
{
	std::vector<int> indicesOfType;
	for (int i = 0; i < m_Frames.m_Frames[frameIndex].size(); i++)
	{
		if (atomLabel == m_Frames.m_Frames[frameIndex][i].m_Symbol)
			indicesOfType.push_back(i);
	}
	return indicesOfType;
}