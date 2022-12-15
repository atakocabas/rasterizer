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

class Scene {
public:
    Color backgroundColor;
    bool cullingEnabled;
    enum REGION {
        FIRST,
        SECOND,
        THIRD,
        FOURTH
    };

    vector<vector<Color> > image;
    vector<Camera *> cameras;
    vector<Vec3 *> vertices;
    vector<Color *> colorsOfVertices;
    vector<Scaling *> scalings;
    vector<Rotation *> rotations;
    vector<Translation *> translations;
    vector<Mesh *> meshes;

    Scene(const char *xmlPath);

    void initializeImage(Camera *camera);

    void forwardRenderingPipeline(Camera *camera);

    int makeBetweenZeroAnd255(double value);

    void writeImageToPPMFile(Camera *camera);

    void convertPPMToPNG(string ppmFileName, int osType);

    Matrix4 ModelingTransformation(Mesh *mesh);

    void triangleRasterization(int i, int j, const vector<vector<Vec3>> &allNewVertexWithVp,Camera *cam );

    void lineRasterization(int i, int j, int id, Camera *cam, vector<vector<Vec3>> vpvertices);

    void interpolate(int x, int y, const Vec3 &a, const Vec3 &b, const Color &color_a, const Color &color_b);

    void drawLine(const Vec3 &smallerVertex, const Vec3 &biggerVertex, REGION region, const Color &color_a,
                  const Color &color_b);

    void writeToImage(const Color &c, int , int );

};

#endif
