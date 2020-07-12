#pragma once

#include <string>
#include <vector>
#include "Vec.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

class Atom
{
public:
	Atom(std::string symbol, glm::vec3 position);

	const std::string m_Symbol;
	glm::vec3 m_Position;
};


class Frames
{
public:
	Frames(const std::string& xyz_file);

	double getFrameVolume(int iFrame);
	std::vector<std::vector<Atom>> m_Frames; // atoms for each frame
	std::vector<glm::vec3> m_CellParameters; // cell parameters for each frame
private:
	void getAtomsAndCell(); // reads xyz file to get each atom and cell parameters
	
	std::string m_XYZFile;
	void getHeaderInfo(std::ifstream& infile, std::string& line);
	bool isInteger(std::string& line);
};