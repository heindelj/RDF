#include "Frames.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem> // exists

Atom::Atom(std::string symbol, Vec<double> position) :
	m_Symbol(symbol), m_Position(position){ }

// constructor which reads xyz file with cell parameters
Frames::Frames(const std::string & xyz_file)
	: m_XYZFile(xyz_file)
{
	getAtomsAndCell();
}


void Frames::getHeaderInfo(std::ifstream& infile, std::string& line)
{
	//std::getline(infile, line); // number of atoms line
	std::istringstream iss(line);

	// these are the main ways I expect the header of the input file might look
	int numAtoms;
	double a, b, c;
	double alpha, beta, gamma, garbage;
	if (iss >> numAtoms >> a >> b >> c)
	{
		m_CellParameters.push_back(Vec<double>(a, b, c));
		std::getline(infile, line); // get comment line of header
		return;
	}
	else if (iss >> numAtoms >> a >> b >> c >> alpha >> beta >> gamma)
	{
		m_CellParameters.push_back(Vec<double>(a, b, c));
		std::getline(infile, line);
		return;
	}
	else if (iss >> numAtoms >> a >> b >> c >> alpha >> beta >> gamma >> garbage)
	{
		m_CellParameters.push_back(Vec<double>(a, b, c));
		std::getline(infile, line);
		return;
	}
	else if (iss >> numAtoms)
	{
		std::getline(infile, line);
		std::istringstream iss(line);
		if (iss >> a >> b >> c)
		{
			m_CellParameters.push_back(Vec<double>(a, b, c));
			return;
		}
	}
	return;
};

void Frames::getAtomsAndCell()
{
	if (!std::filesystem::exists(m_XYZFile))
	{
		throw "Could not find geometry file (.xyz)!";
	}
	//open the file
	std::ifstream infile(m_XYZFile);
	std::string line = "";
	std::string atom_label = "";
	std::vector<Atom> Atoms;

	// temporary variables
	int numAtoms;
	double x, y, z;

	// get initial header
	std::getline(infile, line); 
	getHeaderInfo(infile, line);
	//read all of the coordinates and atom labels for all frames
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		iss >> atom_label >> x >> y >> z;
		if (!isInteger(atom_label))
			Atoms.push_back(Atom(atom_label, Vec<double>(x, y, z)));
		else
		{
			m_Frames.push_back(Atoms);
			Atoms.clear();
			getHeaderInfo(infile, line);
		}
	}
	m_Frames.push_back(Atoms);
}

bool Frames::isInteger(std::string& line)
{
	std::string::const_iterator it = line.begin();
	while (it != line.end() && std::isdigit(*it)) ++it;
	return !line.empty() && it == line.end();
}

double Frames::getFrameVolume(int iFrame)
{
	return m_CellParameters[iFrame].elementProduct();
}

