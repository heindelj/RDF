#include "RDF.h"
#include <cstdlib>     // abs
#define _USE_MATH_DEFINES
#include <math.h>      // sqrt, pi
#include <glm/gtx/component_wise.inl> // compMin

RDF::RDF(Frames frames, int nbins, double maxCutoff) :
	m_Frames(frames), m_NumBins(nbins), m_MaxCutoff(maxCutoff)
{
	std::vector<double> histogram(m_NumBins);
	m_Histogram = histogram;

	// fill in the frame_by_frame_histrogram
	for (int i = 0; i < m_Frames.m_Frames.size(); i++)
		m_FrameByFrameHistogram.push_back(m_Histogram);
	numFrames = m_FrameByFrameHistogram.size();

	m_BinWidth = m_MaxCutoff / double(m_NumBins);
}

void RDF::writeRDF()
{
	// zeroth index will never be populated (and is smaller than any physically meaningful distance)
	for (int i = 0; i < m_Histogram.size(); i++)
		std::cout << m_BinWidth * (double(i) + 1.0) << "   " << m_Histogram[i] << std::endl;
}

void RDF::computeRDF(std::string soluteLabel, std::string solventLabel)
{
	// assumes number of atom type does not change between frames
	numSolventAtoms = getIndicesOfAtomTypeFromFrame(0, solventLabel).size();

	int numThreadedBlocks = int(m_Frames.m_Frames.size() / numThreads);
	for (int iBlock = 0; iBlock < numThreadedBlocks; iBlock++)
	{
		for (int iThread = 0; iThread < numThreads; iThread++)
			threads.emplace_back([=]() {computeSingleFrameRDF(iBlock * numThreads + iThread, soluteLabel, solventLabel, m_Frames.getFrameVolume(iBlock * numThreads + iThread)); });

		for (auto& thread : threads)
			thread.join();

		threads.clear();
	}

	// finish up the ones which didn't fit in a threaded block
	for (int frameIndex = numThreadedBlocks * numThreads; frameIndex < m_Frames.m_Frames.size(); frameIndex++)
		computeSingleFrameRDF(frameIndex, soluteLabel, solventLabel, m_Frames.getFrameVolume(frameIndex));

	// compute average over all frames
	for (int iFrame = 0; iFrame < numFrames; iFrame++)
	{
		for (int iBin = 0; iBin < m_Histogram.size(); iBin++)
			m_Histogram[iBin] += m_FrameByFrameHistogram[iFrame][iBin];
	}


	for (int iBin = 0; iBin < m_Histogram.size(); iBin++)
		m_Histogram[iBin] /= (double)numFrames; // /binwidth
}

void RDF::computeSingleFrameRDF(int frameIndex, std::string soluteLabel, std::string solventLabel, double volume)
{
	if ((frameIndex % 10000) == 0)
		std::cout << "Computing Frame Number: " << frameIndex << std::endl;
	std::vector<int> soluteIndices = getIndicesOfAtomTypeFromFrame(frameIndex, soluteLabel);

	if (soluteLabel == solventLabel)
	{
		for (int iSolute = 0; iSolute < soluteIndices.size() - 1; iSolute++)
		{
			for (int iSolvent = iSolute + 1; iSolvent < soluteIndices.size(); iSolvent++)
			{
				double r = imageDistance(m_Frames.m_Frames[frameIndex][soluteIndices[iSolute]].m_Position, m_Frames.m_Frames[frameIndex][soluteIndices[iSolvent]].m_Position, frameIndex);
				if (r < m_MaxCutoff && r > 0.0)
				{
					int binIndex = int(r / m_BinWidth);
					m_FrameByFrameHistogram[frameIndex][binIndex] += 2.0;
				}
			}
		}
	}
	else
	{
		std::vector<int> solventIndices = getIndicesOfAtomTypeFromFrame(frameIndex, solventLabel);
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
	for (int binIndex = 0; binIndex < m_NumBins; binIndex++)
	{
		m_FrameByFrameHistogram[frameIndex][binIndex] /= ((numSolventAtoms / volume) * 4 * M_PI
			* ((double(binIndex) + 1.0) * m_BinWidth) * ((double(binIndex) + 1.0) * m_BinWidth) * m_BinWidth);
	}
}

double RDF::imageDistance(glm::vec3 vec1, glm::vec3 vec2, int frameIndex)
{
	//glm::vec3 distanceVector = vec1 - vec2;
	// branchless calculation of components
	//double dx = (abs(distanceVector[0])) * (abs(distanceVector[0]) <= (m_Frames.m_CellParameters[frameIndex][0] / 2.0))
	//	+ (m_Frames.m_CellParameters[frameIndex][0] - abs(distanceVector[0])) * (abs(distanceVector[0]) > (m_Frames.m_CellParameters[frameIndex][0] / 2.0));
	//double dy = (abs(distanceVector[1])) * (abs(distanceVector[1]) <= (m_Frames.m_CellParameters[frameIndex][1] / 2.0))
	//	+ (m_Frames.m_CellParameters[frameIndex][1] - abs(distanceVector[1])) * (abs(distanceVector[1]) > (m_Frames.m_CellParameters[frameIndex][1] / 2.0));
	//double dz = (abs(distanceVector[2])) * (abs(distanceVector[2]) <= (m_Frames.m_CellParameters[frameIndex][2] / 2.0))
	//	+ (m_Frames.m_CellParameters[frameIndex][2] - abs(distanceVector[2])) * (abs(distanceVector[2]) > (m_Frames.m_CellParameters[frameIndex][2] / 2.0));
	//return glm::length(glm::vec3(dx, dy, dz));

	return glm::length(vec1 - vec2);
}

double RDF::smallestAbsoluteComponent(glm::vec3 v)
{
	return glm::compMin(glm::vec3(abs(v[0]), abs(v[1]), abs(v[2])));
}

std::vector<int> RDF::getIndicesOfAtomTypeFromFrame(int frameIndex, std::string atomLabel)
{
	std::vector<int> indicesOfType;
	indicesOfType.reserve(m_Frames.m_Frames[frameIndex].size() / 2);
	for (int i = 0; i < m_Frames.m_Frames[frameIndex].size(); i++)
	{
		if (atomLabel == m_Frames.m_Frames[frameIndex][i].m_Symbol)
			indicesOfType.push_back(i);
	}
	return indicesOfType;
}