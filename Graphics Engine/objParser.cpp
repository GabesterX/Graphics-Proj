#include <fstream>
#include <sstream>
#include "GeometryStructs.h"
#include "mtlParser.cpp"

// NOTE - NEXT THING WAS PUSHING TRIANGLES TO SUBMESH TRIGROUP AND THEN PUSHING SUB MESH TO SUB MESH GROUP UPON INSTANTIATION OF NEW SUBMESH OR REACHING END OF FILE


class objParser {
public:

	void parse(const char* filename, std::vector<vec3d> &vBuffer, std::vector<vec2d> &tBuffer, std::vector<subMesh> &meshList) {

		std::ifstream file(filename);
		if (!file) {

			return;
		}
		parse(file, vBuffer, tBuffer, meshList);
	}

	void parse(std::istream& file, std::vector<vec3d> &vBuffer, std::vector<vec2d> &tBuffer, std::vector<subMesh> &meshList) {
		std::string mtlPath;
		std::string line;
		subMesh currMesh;

		while (std::getline(file, line)) {
			std::stringstream ss(line);

			//reads until white space
			ss.unsetf(std::ios_base::skipws);
			ss >> std::ws;

			// If the line is empty
			if (ss.eof()) {
				continue;
			}

			if (ss.peek() == '#') {
				// Comment line in file
				continue;
			}

			std::string keyword;
			ss >> keyword;

			//In the case of a vertex coordinate
			if (keyword == "v") {
				vec3d v;
				ss >> std::ws >> v.x >> std::ws >> v.y >> std::ws >> v.z >> std::ws;
				vBuffer.push_back(v);
				
			}
			//In the case of a texture coordinate
			else if (keyword == "vt") {
				vec2d t;
				ss >> std::ws >> t.u >> std::ws >> t.v >> std::ws;
				t.u = 1.0f - t.u;
				t.v = 1.0f - t.v;
				tBuffer.push_back(t);
			}

			//In the case of a face
			else if (keyword == "f") {
				triangle currFace;
				std::vector<std::string> parts;
				std::vector<int> vertIndices;
				std::vector<int> texIndices;

				while (!ss.eof()) {
					//Read each section is as a part of a list
					int pointcount = 0;
					std::string part;
					ss >> std::ws >> part >> std::ws;
					parts.push_back(part);

					for (unsigned int i = 0; i < parts.size(); i++) {
						int vert;
						int tex;

						std::stringstream partSS(parts[i]);
						partSS >> vert;

						if (!partSS.eof() && (partSS.peek() == '/')) {
							char slash1;
							partSS >> slash1;

							partSS >> tex;
							currFace.texs[pointcount] = tex;
							currFace.texs[pointcount] = currFace.texs[pointcount] - 1;
						}

						currFace.verts[pointcount++] = vert;
						currFace.verts[pointcount - 1] = currFace.verts[pointcount - 1] - 1;
					}
				}
				currMesh.triGroup.push_back(currFace);
			}
			else if (keyword == "o") {
				if (currMesh.name.empty()) {
					ss >> std::ws >> currMesh.name;
				}
				else {
					meshList.push_back(currMesh);
					currMesh = {};
					ss >> std::ws >> currMesh.name;
				}
			}
			else if (keyword == "mtllib") {
				mtlPath = line.substr(7);
			}
			else if (keyword == "usemtl") {
				std::string material;
				ss >> std::ws >> material;
				getTexturePath(material, mtlPath, currMesh);
			}

		}
		meshList.push_back(currMesh);
	}

	void getTexturePath(std::string material, std::string libPath, subMesh &currMesh) {
		std::ifstream file(libPath);
		std::string line;
		std::string texturePath;
		bool foundMat = false;
		while(std::getline(file, line)) {
			std::stringstream ss(line);
			std::string keyword;
			ss.unsetf(std::ios_base::skipws);

			ss >> std::ws;

			ss >> std::ws >> keyword;
			std::string currMat;
			if (foundMat) {
				std::string matName;
				ss >> std::ws >> matName;
				if (keyword.length() >= 3) {
					if (keyword.substr(0, 3) == "map") {
						texturePath = line.substr(6 + (keyword.length() - 5));
						
						currMesh.subTex = new olc::Sprite(texturePath);
						currMesh.hasTexture = true;
						return;
					}
				}
				else if (matName == "") {
					currMesh.hasTexture = false;
					return;
				}
			}

			if (keyword == "newmtl") {
				currMat = line.substr(7);
				if (currMat == material) {
					foundMat = true;
				}
			}
		}
		currMesh.hasTexture = false;
	}
};