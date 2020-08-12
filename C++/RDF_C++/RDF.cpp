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
	binRadii.reserve(m_NumBins);
	for (int i = 0; i <= m_NumBins; i++)
		binRadii.push_back(m_BinWidth * (double(i) + 1.0));

	// fold the images back into the PBC box so distances can be calculated normally.
	//foldImagesIntoBox();
}

void RDF::writeRDF()
{
	// should report the center of the bin to get correct average pair distance
	for (int i = 0; i < m_NumBins; i++)
		std::cout << binRadii[i] - m_BinWidth / 2.0 << "   " << m_Histogram[i] << std::endl;
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
		for (int iBin = 0; iBin < m_NumBins; iBin++)
			m_Histogram[iBin] += m_FrameByFrameHistogram[iFrame][iBin] / (double)numFrames;
	}
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
		m_FrameByFrameHistogram[frameIndex][binIndex] /= (4 * M_PI * ((double)numSolventAtoms * (double)(numSolventAtoms - 1)) / volume
			* (binRadii[binIndex] * binRadii[binIndex]) * m_BinWidth);
	}
}

double RDF::imageDistance(glm::vec3 vec1, glm::vec3 vec2, int frameIndex)
{
	
	glm::vec3 distanceVector = vec1 - vec2;
	const glm::vec3& cellParams = m_Frames.m_CellParameters[frameIndex];
	const glm::vec3 testDistance = cellParams * 0.5f;
	// branchless calculation of components
	
	float dx = (distanceVector[0] + cellParams[0] * int(glm::abs((distanceVector[0] - testDistance[0]) / cellParams[0]))) * (distanceVector[0] < (-testDistance[0]))
		+ (distanceVector[0] - cellParams[0] * int(glm::abs((distanceVector[0] + testDistance[0]) / cellParams[0]))) * (distanceVector[0] > testDistance[0])
		+ distanceVector[0] * (glm::abs(distanceVector[0]) <= testDistance[0]);
	float dy = (distanceVector[1] + cellParams[1] * int(glm::abs((distanceVector[1] - testDistance[1]) / cellParams[1]))) * (distanceVector[1] < (-testDistance[1]))
		+ (distanceVector[1] - cellParams[1] * int(glm::abs((distanceVector[1] + testDistance[1]) / cellParams[1]))) * (distanceVector[1] > testDistance[1])
		+ distanceVector[1] * (glm::abs(distanceVector[1]) <= testDistance[1]);
	float dz = (distanceVector[2] + cellParams[2] * int(glm::abs((distanceVector[2] - testDistance[2]) / cellParams[2]))) * (distanceVector[2] < (-testDistance[2]))
		+ (distanceVector[2] - cellParams[2] * int(glm::abs((distanceVector[2] + testDistance[2]) / cellParams[2]))) * (distanceVector[2] > testDistance[2])
		+ distanceVector[2] * (glm::abs(distanceVector[2]) <= testDistance[2]);
	return glm::length(glm::vec3(dx, dy, dz));
	
	//return glm::length(vec1 - vec2);
}

void RDF::foldImagesIntoBox()
{
	int N = 69;
	std::cout << m_Frames.m_Frames[0][N].m_Position[0] << " " << m_Frames.m_Frames[0][N].m_Position[1] << " " << m_Frames.m_Frames[0][N].m_Position[2] << std::endl;
	for (int iFrame = 0; iFrame < m_Frames.m_Frames.size(); iFrame++)
	{
		const glm::vec3& cellParams = m_Frames.m_CellParameters[iFrame];
		const glm::vec3 testDistance = cellParams * 0.5f;
		for (int iAtom = 0; iAtom < m_Frames.m_Frames[iFrame].size(); iAtom++)
		{
			glm::vec3& pos = m_Frames.m_Frames[iFrame][iAtom].m_Position;
			for (int w = 0; w < 3; w++)
			{
				m_Frames.m_Frames[iFrame][iAtom].m_Position[w] += cellParams[w] * int(glm::abs(nearbyint(pos[w] / cellParams[w]))) * (pos[w] < -testDistance[w])
					- cellParams[w] * int(glm::abs(nearbyint(pos[w] / cellParams[w]))) * (pos[w] > testDistance[w]);
			}
		}
	}
	std::cout << m_Frames.m_Frames[0][N].m_Position[0] << " " << m_Frames.m_Frames[0][N].m_Position[1] << " " << m_Frames.m_Frames[0][N].m_Position[2] << std::endl;
}

double RDF::smallestAbsoluteComponent(glm::vec3 v)
{
	return glm::compMin(glm::vec3(abs(v[0]), abs(v[1]), abs(v[2])));
}

std::vector<int> RDF::getIndicesOfAtomTypeFromFrame(int frameIndex, std::string& atomLabel)
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
