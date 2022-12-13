#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Matrix4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;
    enum REGION {
        FIRST,
        SECOND,
        THIRD,
        FOURTH
    };

    vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Mesh* > meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);
    Matrix4 ModelingTransformation(Mesh* mesh);
	Matrix4 ViewingTransformation(Mesh *mesh);
    void rasterization(int, int, vector<vector<Vec3>>);
	void midPoint(int i, int j, int id, Camera* cam, vector< vector<Vec3> > vpvertices);
	void draw(int x, int y, Vec3 a, Vec3 b);
	void drawLine( Vec3 a, Vec3 b,REGION region);
};

#endif
