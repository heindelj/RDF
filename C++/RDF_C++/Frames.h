#pragma once

#include <string>
#include <vector>
#include "Vec.h"

class Atom
{
public:
	Atom(std::string symbol, Vec<double> position);

	const std::string m_Symbol;
	Vec<double> m_Position;

	// store possible atom pairs here so that you can avoid index-ordering dependence
	std::vector<int> m_PairList;
};


class Frames
{
public:
	Frames(const std::string& xyz_file);

	double getFrameVolume(int iFrame);
	std::vector<std::vector<Atom>> m_Frames; // atoms for each frame
	std::vector<Vec<double>> m_CellParameters; // cell parameters for each frame
private:
	void getAtomsAndCell(); // reads xyz file to get each atom and cell parameters
	
	std::string m_XYZFile;
	void getHeaderInfo(std::ifstream& infile, std::string& line);
	bool isInteger(std::string& line);
};