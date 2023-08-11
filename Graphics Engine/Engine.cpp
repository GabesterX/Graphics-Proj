#define OLC_PGE_APPLICATION
#define _SILENCE_CXX17_STRSTREAM_DEPRECATION_WARNING
#define _CRT_SECURE_NO_WARNINGS
//#include "olcPixelGameEngine.h"
//#include "GeometryStructs.h"
#include "objParser.cpp"
#include <fstream>
#include <strstream>
#include <algorithm>




struct mesh {
	std::string mtllibPath;
	std::vector<triangle> tris;
	std::vector<subMesh> meshList;
	/*
	bool LoadFromObjectFile(std::string sFilename, std::vector<vec3d>& out_vBuffer, std::vector<vec2d>& out_tBuffer, bool bHasTexture = true)
	{
		std::ifstream f(sFilename);
		if (!f.is_open())
			return false;

		// Local cache of verts
		std::vector<vec3d> verts;
		std::vector<vec2d> texs;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				if (line[1] == 't')
				{
					triangle currTri;
					vec2d v;
					s >> junk >> junk >> v.u >> v.v;
					// A little hack for the spyro texture
					v.u = 1.0f - v.u;
					v.v = 1.0f - v.v;
					out_tBuffer.push_back(v);
				}
				else
				{
					vec3d v;
					s >> junk >> v.x >> v.y >> v.z;
					out_vBuffer.push_back(v);
				}
			}

			if (!bHasTexture)
			{
				if (line[0] == 'f')
				{

					triangle currTri;

					int f[3];
					s >> junk >> f[0] >> f[1] >> f[2];
					currTri.verts[0] = f[0];
					currTri.verts[1] = f[1];
					currTri.verts[2] = f[2];
					tris.push_back(currTri);
				}
			}
			else
			{
				if (line[0] == 'f')
				{

					triangle currTri;

					s >> junk;


					std::string tokens[6];
					int nTokenCount = -1;


					while (!s.eof())
					{
						char c = s.get();
						if (c == ' ' || c == '/')
							nTokenCount++;
						else
							tokens[nTokenCount].append(1, c);
					}

					tokens[nTokenCount].pop_back();

					currTri.verts[0] = stoi(tokens[0]) - 1;
					currTri.texs[0] = stoi(tokens[1]) - 1;
					currTri.verts[1] = stoi(tokens[2]) - 1;
					currTri.texs[1] = stoi(tokens[3]) - 1;
					currTri.verts[2] = stoi(tokens[4]) - 1;
					currTri.texs[2] = stoi(tokens[5]) - 1;

					tris.push_back(currTri);

				}

			}
		}
		return true;
	}
	*/
};






// Override base class with your custom functionality
class epqGraphics : public olc::PixelGameEngine
{
public:
	epqGraphics()
	{
		// Name your application
		sAppName = "EPQ Graphics";
	}

private:
	std::vector<subMesh> meshList;
	//mesh meshObj;
	mat4d matProj;
	vec3d vCamera;
	vec3d vLookDir;
	float fTheta;
	float fYaw;
	olc::Sprite *sprTex1;
	
	std::vector<vec3d> vertBuffer;
	std::vector<vec2d> texBuffer;

	vec3d Matrix_MultiplyVector(mat4d& m, vec3d& i) {
		vec3d v;
		v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
		v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
		v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
		v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
		return v;
	}

	mat4d Matrix_MakeIdentity() {
		mat4d matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4d Matrix_MakeRotationX(float fAngleRad) {
		mat4d matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[1][2] = sinf(fAngleRad);
		matrix.m[2][1] = -sinf(fAngleRad);
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4d Matrix_MakeRotationY(float fAngleRad) {
		mat4d matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][2] = sinf(fAngleRad);
		matrix.m[2][0] = -sinf(fAngleRad);
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4d Matrix_MakeRotationZ(float fAngleRad) {
		mat4d matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][1] = sinf(fAngleRad);
		matrix.m[1][0] = -sinf(fAngleRad);
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	mat4d Matrix_MakeTranslation(float x, float y, float z) {
		mat4d matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		matrix.m[3][0] = x;
		matrix.m[3][1] = y;
		matrix.m[3][2] = z;
		return matrix;
	}

	mat4d Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar) {
		float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
		mat4d matrix;
		matrix.m[0][0] = fAspectRatio * fFovRad;
		matrix.m[1][1] = fFovRad;
		matrix.m[2][2] = fFar / (fFar - fNear);
		matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
		matrix.m[2][3] = 1.0f;
		matrix.m[3][3] = 0.0f;
		return matrix;
	}

	mat4d Matrix_MultiplyMatrix(mat4d& m1, mat4d& m2) {
		mat4d matrix;
		for (int c = 0; c < 4; c++) {
			for (int r = 0; r < 4; r++) {
				matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
			}
		}
		return matrix;
	}

	vec3d Vector_Add(vec3d& v1, vec3d& v2) {
		return  { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
	}

	vec3d Vector_Sub(vec3d& v1, vec3d& v2) {
		return  { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
	}

	vec3d Vector_Mul(vec3d& v1, float k)
	{
		return { v1.x * k, v1.y * k, v1.z * k };
	}


	vec3d Vector_Div(vec3d& v1, float k) {
		return { v1.x / k, v1.y / k, v1.z / k };
	}

	float Vector_DotProduct(vec3d& v1, vec3d& v2) {
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	mat4d Matrix_PointAt(vec3d& pos, vec3d& target, vec3d& up) {
		// Calc new forward direction
		vec3d newForward = Vector_Sub(target, pos);
		newForward = Vector_Normalise(newForward);

		//Calculate new UP direction
		vec3d a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));
		vec3d newUp = Vector_Sub(up, a);
		newUp = Vector_Normalise(newUp);

		// New Right Direction is Cross product
		vec3d newRight = Vector_CrossProduct(newUp, newForward);

		// Dimensioning and translation Matrix
		mat4d matrix;
		matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
		matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
		return matrix;

	}

	mat4d Matrix_QuickInverse(mat4d& m) // Only for Rotation/Translation Matrices
	{
		mat4d matrix;
		matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
		matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
		matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
		matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}

	float Vector_Length(vec3d& v) {
		return sqrtf(Vector_DotProduct(v, v));
	}

	vec3d Vector_Normalise(vec3d& v) {
		float l = Vector_Length(v);
		return { v.x / l, v.y / l, v.z / l };
	}


	vec3d Vector_CrossProduct(vec3d& v1, vec3d& v2) {
		vec3d v;
		v.x = v1.y * v2.z - v1.z * v2.y;
		v.y = v1.z * v2.x - v1.x * v2.z;
		v.z = v1.x * v2.y - v1.y * v2.x;
		return v;
	}

	vec3d Vector_IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd, float &t) {
		plane_n = Vector_Normalise(plane_n);
		float plane_d = -Vector_DotProduct(plane_n, plane_p);
		float ad = Vector_DotProduct(lineStart, plane_n);
		float bd = Vector_DotProduct(lineEnd, plane_n);
		t = (-plane_d - ad) / (bd - ad);
		vec3d lineStartToEnd = Vector_Sub(lineEnd, lineStart);
		vec3d lineToIntersect = Vector_Mul(lineStartToEnd, t);
		return Vector_Add(lineStart, lineToIntersect);
	}

	int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, primitive& in_tri, primitive& out_tri1, primitive& out_tri2)
	{
		// Make sure plane normal is indeed normal
		plane_n = Vector_Normalise(plane_n);

		// Return signed shortest distance from point to plane, plane normal must be normalised
		auto dist = [&](vec3d& p)
		{
			vec3d n = Vector_Normalise(p);
			return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
		};

		// Create two temporary storage arrays to classify points either side of plane
		// If distance sign is positive, point lies on "inside" of plane
		vec3d* inside_points[3];  int nInsidePointCount = 0;
		vec3d* outside_points[3]; int nOutsidePointCount = 0;

		vec2d* inside_tex[3]; int nInsideTexCount = 0;
		vec2d* outside_tex[3]; int nOutsideTexCount = 0;

		// Get signed distance of each point in triangle to plane
		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; inside_tex[nInsideTexCount++] = &in_tri.t[0]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; outside_tex[nOutsideTexCount++] = &in_tri.t[0]; }
		if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; inside_tex[nInsideTexCount++] = &in_tri.t[1]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; outside_tex[nOutsideTexCount++] = &in_tri.t[1]; }
		if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; inside_tex[nInsideTexCount++] = &in_tri.t[2]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; outside_tex[nOutsideTexCount++] = &in_tri.t[2]; }

		// Now classify triangle points, and break the input triangle into
		// smaller output triangles if required. There are four possible
		// outcomes...

		if (nInsidePointCount == 0)
		{
			// All points lie on the outside of plane, so clip whole triangle
			// It ceases to exist

			return 0; // No returned triangles are valid
		}

		if (nInsidePointCount == 3)
		{
			// All points lie on the inside of plane, so do nothing
			// and allow the triangle to simply pass through
			out_tri1 = in_tri;

			return 1; // Just the one returned original triangle is valid
		}

		if (nInsidePointCount == 1 && nOutsidePointCount == 2)
		{
			// Triangle should be clipped. As two points lie outside
			// the plane, the triangle simply becomes a smaller triangle

			// Copy appearance info to new triangle


			// The inside point is valid, so keep that...
			out_tri1.p[0] = *inside_points[0];
			out_tri1.t[0] = *inside_tex[0];
			// but the two new points are at the locations where the
			// original sides of the triangle (lines) intersect with the plane
			
			// t for intersection point when triangle is clipped, required to calculate new texture coordinates.
			float t;
			out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0],t);
			out_tri1.t[1].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[1].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[1].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;


			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1],t);
			out_tri1.t[2].u = t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[2].w = t * (outside_tex[1]->w - inside_tex[0]->w) + inside_tex[0]->w;

			return 1; // Return the newly formed single triangle
		}

		if (nInsidePointCount == 2 && nOutsidePointCount == 1)
		{
			// Triangle should be clipped. As two points lie inside the plane,
			// the clipped triangle becomes a "quad". Fortunately, we can
			// represent a quad with two new triangles

			// Copy appearance info to new triangles
			// The first triangle consists of the two inside points and a new
			// point determined by the location where one side of the triangle
			// intersects with the plane
			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.t[0] = *inside_tex[0];
			out_tri1.t[1] = *inside_tex[1];


			float t;
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0],t );
			out_tri1.t[2].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[2].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			// The second triangle is composed of one of he inside points, a
			// new point determined by the intersection of the other side of the
			// triangle and the plane, and the newly created point above
			out_tri2.p[0] = *inside_points[1];
			out_tri2.t[0] = *inside_tex[1];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.t[1] = out_tri1.t[2];
			out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0],t);
			out_tri2.t[2].u = t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
			out_tri2.t[2].v = t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;
			out_tri2.t[2].w = t * (outside_tex[0]->w - inside_tex[1]->w) + inside_tex[1]->w;

			return 2; // Return two newly formed triangles which form a quad
		}
	}


	olc::Pixel GetColour(float lum) {

		int nValue = (int)(std::max(lum, 0.20f) * 255.0f);
		return olc::Pixel(nValue, nValue, nValue);

	}

	float* pDepthBuffer = nullptr;

public:
	bool OnUserCreate() override
	{
		// Called once at the start, so create things here
		
		

		pDepthBuffer = new float[ScreenWidth() * ScreenHeight()];
		
		

		objParser parser;
		parser.parse("dreamland.obj", vertBuffer, texBuffer, meshList);

		

		//meshObj.LoadFromObjectFile("SpyroLevelOne.obj", vertBuffer, texBuffer);


		sprTex1 = new olc::Sprite("High2.png");

		

		//Projection Matrix
		matProj = Matrix_MakeProjection(90.0f, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.f);

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{

		if (GetKey(olc::Key::SPACE).bHeld) {
			vCamera.y += 8.0f * fElapsedTime;
		}

		if (GetKey(olc::Key::CTRL).bHeld) {
			vCamera.y -= 8.0f * fElapsedTime;
		}

		if (GetKey(olc::Key::LEFT).bHeld) {
			vCamera.x += 8.0f * fElapsedTime;	// Travel Along X-Axis
		}

		if (GetKey(olc::Key::RIGHT).bHeld) {
			vCamera.x -= 8.0f * fElapsedTime;	// Travel Along X-Axis
		}

		vec3d vForward = Vector_Mul(vLookDir, 16.0f * fElapsedTime);

		// Standard FPS Control scheme, but turn instead of strafe
		if (GetKey(olc::Key::W).bHeld) {
			vCamera = Vector_Add(vCamera, vForward);
		}
		if (GetKey(olc::Key::S).bHeld) {
			vCamera = Vector_Sub(vCamera, vForward);
		}
		if (GetKey(olc::Key::A).bHeld) {
			fYaw -= 2.0f * fElapsedTime;
		}
		if (GetKey(olc::Key::D).bHeld) {
			fYaw += 2.0f * fElapsedTime;
		}





		mat4d matRotZ, matRotX;
		matRotZ = Matrix_MakeRotationZ(fTheta * 0.5f);
		matRotX = Matrix_MakeRotationX(fTheta);
		mat4d matTrans;
		matTrans = Matrix_MakeTranslation(0.0f, 0.0f, 16.0f);
		mat4d matWorld;
		matWorld = Matrix_MakeIdentity();
		matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);
		matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);

		vec3d vUp = { 0,1,0 };
		vec3d vTarget = { 0,0,1 };
		mat4d matCameraRot = Matrix_MakeRotationY(fYaw);
		vLookDir = Matrix_MultiplyVector(matCameraRot, vTarget);
		vTarget = Vector_Add(vCamera, vLookDir);
		mat4d matCamera = Matrix_PointAt(vCamera, vTarget, vUp);
		mat4d matView = Matrix_QuickInverse(matCamera);

		std::vector<primitive> vecTrianglesToRaster;

		for (auto& currMesh : meshList) {
			
			for (auto& tri : currMesh.triGroup) {


				primitive currTri;


				// World transform

				currTri.p[0] = Matrix_MultiplyVector(matWorld, vertBuffer[tri.verts[0]]);
				currTri.p[1] = Matrix_MultiplyVector(matWorld, vertBuffer[tri.verts[1]]);
				currTri.p[2] = Matrix_MultiplyVector(matWorld, vertBuffer[tri.verts[2]]);

				currTri.t[0] = texBuffer[tri.texs[0]];
				currTri.t[1] = texBuffer[tri.texs[1]];
				currTri.t[2] = texBuffer[tri.texs[2]];

				//Use Cross-Product to get surface normal
				vec3d normal, L1, L2;

				//Lines either side of the triangle
				L1 = Vector_Sub(currTri.p[1], currTri.p[0]);
				L2 = Vector_Sub(currTri.p[2], currTri.p[0]);

				// Take cross product of lines to get normal to triangle surface
				normal = Vector_CrossProduct(L1, L2);

				// Normalise the line

				normal = Vector_Normalise(normal);


				//Get Ray from Triangle to Camera

				vec3d vCameraRay = Vector_Sub(currTri.p[0], vCamera);

				if (Vector_DotProduct(normal, vCameraRay) < 0.0f) {

					/*Illumination
					vec3d light_direction = { 0.0f, 1.0f, -1.0f };
					light_direction = Vector_Normalise(light_direction);

					// How "Aligned" are light direction and triangle surface normal?

					float dp = std::max(0.1f, Vector_DotProduct(light_direction, normal));


					currTri.col = GetColour(dp);
					*/

					// Convert World Space to View Space
					currTri.p[0] = Matrix_MultiplyVector(matView, currTri.p[0]);
					currTri.p[1] = Matrix_MultiplyVector(matView, currTri.p[1]);
					currTri.p[2] = Matrix_MultiplyVector(matView, currTri.p[2]);

					currTri.t[0] = currTri.t[0];
					currTri.t[1] = currTri.t[1];
					currTri.t[2] = currTri.t[2];

					int nClippedTriangles = 0;
					primitive clipped[2];
					nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, currTri, clipped[0], clipped[1]);

					for (int n = 0; n < nClippedTriangles; n++) {


						// PROJECT TRIANGLES FROM 3D --> 2D
						currTri.p[0] = Matrix_MultiplyVector(matProj, clipped[n].p[0]);
						currTri.p[1] = Matrix_MultiplyVector(matProj, clipped[n].p[1]);
						currTri.p[2] = Matrix_MultiplyVector(matProj, clipped[n].p[2]);
						currTri.t[0] = clipped[n].t[0];
						currTri.t[1] = clipped[n].t[1];
						currTri.t[2] = clipped[n].t[2];

						currTri.t[0].u = currTri.t[0].u / currTri.p[0].w;
						currTri.t[1].u = currTri.t[1].u / currTri.p[1].w;
						currTri.t[2].u = currTri.t[2].u / currTri.p[2].w;

						currTri.t[0].v = currTri.t[0].v / currTri.p[0].w;
						currTri.t[1].v = currTri.t[1].v / currTri.p[1].w;
						currTri.t[2].v = currTri.t[2].v / currTri.p[2].w;

						currTri.t[0].w = 1.0f / currTri.p[0].w;
						currTri.t[1].w = 1.0f / currTri.p[1].w;
						currTri.t[2].w = 1.0f / currTri.p[2].w;

						//Perspective Divide
						currTri.p[0] = Vector_Div(currTri.p[0], currTri.p[0].w);
						currTri.p[1] = Vector_Div(currTri.p[1], currTri.p[1].w);
						currTri.p[2] = Vector_Div(currTri.p[2], currTri.p[2].w);


						currTri.p[0].x *= -1.0f;
						currTri.p[1].x *= -1.0f;
						currTri.p[2].x *= -1.0f;
						currTri.p[0].y *= -1.0f;
						currTri.p[1].y *= -1.0f;
						currTri.p[2].y *= -1.0f;




						//Offset Verts into visible normalised space
						vec3d vOffsetView = { 1,1,0 };

						currTri.p[0] = Vector_Add(currTri.p[0], vOffsetView);
						currTri.p[1] = Vector_Add(currTri.p[1], vOffsetView);
						currTri.p[2] = Vector_Add(currTri.p[2], vOffsetView);


						currTri.p[0].x *= 0.5 * (float)ScreenWidth();
						currTri.p[0].y *= 0.5 * (float)ScreenHeight();
						currTri.p[1].x *= 0.5 * (float)ScreenWidth();
						currTri.p[1].y *= 0.5 * (float)ScreenHeight();
						currTri.p[2].x *= 0.5 * (float)ScreenWidth();
						currTri.p[2].y *= 0.5 * (float)ScreenHeight();

						//Store Triangle for Sorting (Painters Algorithm)
						vecTrianglesToRaster.push_back(currTri);

						/*
						FillTriangle(currTri.p[0].x, currTri.p[0].y, currTri.p[1].x, currTri.p[1].y, currTri.p[2].x, currTri.p[2].y, currTri.col);
						DrawTriangle(currTri.p[0].x, currTri.p[0].y, currTri.p[1].x, currTri.p[1].y, currTri.p[2].x, currTri.p[2].y, olc::BLACK);
						*/
					}
				}
			}

			/*
			//Sort triangles from back to front
			sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(), [](triangle& t1, triangle& t2) {
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
				return z1 > z2;
				});
			*/



			for (auto& triToRaster : vecTrianglesToRaster) {

				//Clip triangle against 4 edges of the screen
				primitive clipped[2];

				std::list<primitive> listTriangles;

				listTriangles.push_back(triToRaster);

				int nNewTriangles = 1;

				for (int p = 0; p < 4; p++) {
					int nTrisToAdd = 0;
					while (nNewTriangles > 0) {
						primitive test = listTriangles.front();
						listTriangles.pop_front();
						nNewTriangles--;
						switch (p) {
						case 0: nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
						case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
						case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
						case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
						}
						for (int w = 0; w < nTrisToAdd; w++) {
							listTriangles.push_back(clipped[w]);
						}

					}
					nNewTriangles = listTriangles.size();
				}

				for (auto& t : listTriangles) {

					if (currMesh.hasTexture) {
						TexturedTriangle(t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, (1.0f / t.p[0].z), t.t[0].w,
							t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, (1.0f / t.p[1].z), t.t[1].w,
							t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, (1.0f / t.p[2].z), t.t[2].w, currMesh.subTex);
						DrawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, olc::WHITE);
					}
					else {
						DrawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, olc::WHITE);
					}

					

					//FillTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, t.col);
					//DrawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, olc::WHITE);
				}

			}

			

			
		}
		return true;
	}
	void TexturedTriangle(int x1, int y1, float u1, float v1, float z1, float w1,
						  int x2, int y2, float u2, float v2, float z2, float w2,
						  int x3, int y3, float u3, float v3, float z3, float w3,
						  olc::Sprite *tex) {
		if (y2 < y1) {
			std::swap(y1, y2);
			std::swap(x1, x2);
			std::swap(u1, u2);
			std::swap(v1, v2);
			std::swap(z1, z2);
			std::swap(w1, w2);
		}

		if (y3 < y1) {
			std::swap(y1, y3);
			std::swap(x1, x3);
			std::swap(u1, u3);
			std::swap(v1, v3);
			std::swap(z1, z3);
			std::swap(w1, w3);
		}

		if (y3 < y2) {
			std::swap(y2, y3);
			std::swap(x2, x3);
			std::swap(u2, u3);
			std::swap(v2, v3);
			std::swap(z2, z3);
			std::swap(w2, w3);
		}

		int dy1 = y2 - y1;
		int dx1 = x2 - x1;
		float dv1 = v2 - v1;
		float du1 = u2 - u1;
		float dz1 = z2 - z1;
		float dw1 = w2 - w1;


		int dy2 = y3 - y1;
		int dx2 = x3 - x1;
		float dv2 = v3 - v1;
		float du2 = u3 - u1;
		float dz2 = z3 - z1;
		float dw2 = w3 - w1;

		float tex_u, tex_v, tex_z, tex_w;

		float dax_step = 0, dbx_step = 0, du1_step = 0, dv1_step = 0, du2_step = 0, dv2_step = 0, dz1_step = 0, dz2_step = 0, dw1_step = 0, dw2_step = 0;

		if (dy1) dax_step = dx1 / (float)abs(dy1);
		if (dy2) dbx_step = dx2 / (float)abs(dy2);

		if (dy1) du1_step = du1 / (float)abs(dy1);
		if (dy1) dv1_step = dv1 / (float)abs(dy1);
		if (dy1) dz1_step = dz1 / (float)abs(dy1);
		if (dy1) dw1_step = dw1 / (float)abs(dy1);

		if (dy2) du2_step = du2 / (float)abs(dy2);
		if (dy2) dv2_step = dv2 / (float)abs(dy2);
		if (dy2) dz2_step = dz2 / (float)abs(dy2);
		if (dy2) dw2_step = dw2 / (float)abs(dy2);

		if (dy1) {
			for (int i = y1; i <= y2; i++) {
				
				int ax = x1 + (float)(i - y1) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				float tex_su = u1 + (float)(i - y1) * du1_step;
				float tex_sv = v1 + (float)(i - y1) * dv1_step;
				float tex_sz = z1 + (float)(i - y1) * dz1_step;
				float tex_sw = w1 + (float)(i - y1) * dw1_step;

				float tex_eu = u1 + (float)(i - y1) * du2_step;
				float tex_ev = v1 + (float)(i - y1) * dv2_step;
				float tex_ez = z1 + (float)(i - y1) * dz2_step;
				float tex_ew = w1 + (float)(i - y1) * dw2_step;

				if (ax > bx) {
					std::swap(ax, bx);
					std::swap(tex_su, tex_eu);
					std::swap(tex_sv, tex_ev);
					std::swap(tex_sz, tex_ez);
					std::swap(tex_sw, tex_ew);
				}

				tex_u = tex_su;
				tex_v = tex_sv;
				tex_z = tex_sz;
				tex_w = tex_sw;

				float tstep = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				for (int j = ax; j < bx; j++) {

					tex_u = (1.0f - t) * tex_su + t * tex_eu;
					tex_v = (1.0f - t) * tex_sv + t * tex_ev;
					tex_z = (1.0f - t) * tex_sz + t * tex_ez;
					tex_w = (1.0f - t) * tex_sw + t * tex_ew;
					if (tex_z > pDepthBuffer[i * ScreenWidth() + j])
					{
						Draw(j, i, tex->Sample(tex_u / tex_w, tex_v / tex_w));
						pDepthBuffer[i * ScreenWidth() + j] = tex_z;
					}
					t += tstep;
				}
				

			}

		}

		dy1 = y3 - y2;
		dx1 = x3 - x2;
		dv1 = v3 - v2;
		du1 = u3 - u2;
		dz1 = z3 - z2;
		dw1 = w3 - w2;

		if (dy1) dax_step = dx1 / (float)abs(dy1);
		if (dy2) dbx_step = dx2 / (float)abs(dy2);

		du1_step = 0, dv1_step = 0;

		if (dy1) du1_step = du1 / (float)abs(dy1);
		if (dy1) dv1_step = dv1 / (float)abs(dy1);
		if (dy1) dz1_step = dz1 / (float)abs(dy1);
		if (dy1) dw1_step = dw1 / (float)abs(dy1);
		
		if (dy1){
			for (int i = y2; i <= y3; i++) {
				
				int ax = x2 + (float)(i - y2) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				float tex_su = u2 + (float)(i - y2) * du1_step;
				float tex_sv = v2 + (float)(i - y2) * dv1_step;
				float tex_sz = z2 + (float)(i - y2) * dz1_step;
				float tex_sw = w2 + (float)(i - y2) * dw1_step;

				float tex_eu = u1 + (float)(i - y1) * du2_step;
				float tex_ev = v1 + (float)(i - y1) * dv2_step;
				float tex_ez = z1 + (float)(i - y1) * dz2_step;
				float tex_ew = w1 + (float)(i - y1) * dw2_step;


				if (ax > bx) {
					std::swap(ax, bx);
					std::swap(tex_su, tex_eu);
					std::swap(tex_sv, tex_ev);
					std::swap(tex_sz, tex_ez);
					std::swap(tex_sw, tex_ew);
				}

				tex_u = tex_su;
				tex_v = tex_sv;
				tex_z = tex_sz;
				tex_w = tex_sw;


				float tstep = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				for (int j = ax; j < bx; j++) {
					tex_u = (1.0f - t) * tex_su + t * tex_eu;
					tex_v = (1.0f - t) * tex_sv + t * tex_ev;
					tex_z = (1.0f - t) * tex_sz + t * tex_ez;
					tex_w = (1.0f - t) * tex_sw + t * tex_ew;

					if (tex_z > pDepthBuffer[i * ScreenWidth() + j])
					{
						Draw(j, i, tex->Sample(tex_u / tex_w, tex_v / tex_w));
						pDepthBuffer[i * ScreenWidth() + j] = tex_z;
					}
					t += tstep;
				}

			}


		}
	}


};


int main()
{
	epqGraphics demo;
	if (demo.Construct(1366, 768, 4, 4)){
		demo.Start();
	}

	return 0;
}