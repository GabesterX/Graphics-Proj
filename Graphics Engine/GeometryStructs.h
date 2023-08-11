#pragma once
#include <vector>
#include "olcPixelGameEngine.h"

struct vec2d {
	float u = 0;
	float v = 0;
	float w = 1.0f;
};

struct vec3d {
	float x = 0;
	float y = 0;
	float z = 0;
	float w = 1;
};

struct triangle {
	int verts[3];
	int texs[3];
};

struct primitive {
	vec3d p[3];
	vec2d t[3];
};

struct mat4d {
	float m[4][4] = { 0 }; // Row - Column - Initialise to 0
};

struct subMesh {
	bool hasTexture = false;
	std::string name;
	std::string material;
	olc::Sprite* subTex = new olc::Sprite();
	std::vector<triangle> triGroup;
};