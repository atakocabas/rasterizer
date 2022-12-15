#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

void Scene::forwardRenderingPipeline(Camera *camera) {

    vector<vector<Vec3>> allNewVertex;
    vector<vector<Vec3>> allNewVertexWithVp;
    Matrix4 orthPerspectiveProjectionMatrix;
    Matrix4 camViewMatrix = computeViewingTransformationMatrix(camera);
    if (camera->projectionType == 0) { // ORTHOGRAPHIC
        orthPerspectiveProjectionMatrix = calculateOrthographicProjection(camera);
    } else {
        orthPerspectiveProjectionMatrix = calculatePerspectiveProjection(camera);
    }
    Matrix4 viewportTransformationMatrix = calculateViewportMatrix(camera);
    for (auto &mesh: this->meshes) {
        Matrix4 modelingViewMatrix = ModelingTransformation(mesh);
        modelingViewMatrix = multiplyMatrixWithMatrix(orthPerspectiveProjectionMatrix,
                                                      multiplyMatrixWithMatrix(camViewMatrix, modelingViewMatrix));
        vector<Vec3> meshVertex;
        vector<Vec3> meshVertexWithVp;

        for (auto &triangle: mesh->triangles) {
            for (auto &vertex_id: triangle.vertexIds) {
                // TODO: Check Ids if any range error thrown
                Vec3 *new_vertex = this->vertices[vertex_id - 1];
                Vec4 vertex = {
                        new_vertex->x,
                        new_vertex->y,
                        new_vertex->z,
                        1,
                        new_vertex->colorId,
                };

                vertex = multiplyMatrixWithVec4(modelingViewMatrix, vertex);
                vertex.x = vertex.x / vertex.t;
                vertex.y = vertex.y / vertex.t;
                vertex.z = vertex.z / vertex.t;
                vertex.t = 1;

                // ViewPort here?
                Vec3 mshVertex = {
                        vertex.x,
                        vertex.y,
                        vertex.z,
                        vertex.colorId,
                };

                Vec3 mshVertexWithVp = {
                        // vertex ile viewport matrixi çarpacağız.
                        (vertex.x * viewportTransformationMatrix.val[0][0]) +
                        (vertex.y * viewportTransformationMatrix.val[0][1]) +
                        (vertex.z * viewportTransformationMatrix.val[0][2]) +
                        (vertex.t * viewportTransformationMatrix.val[0][3]),

                        ((vertex.x * viewportTransformationMatrix.val[1][0]) +
                         (vertex.y * viewportTransformationMatrix.val[1][1]) +
                         (vertex.z * viewportTransformationMatrix.val[1][2]) +
                         (vertex.t * viewportTransformationMatrix.val[1][3])),

                        ((vertex.x * viewportTransformationMatrix.val[2][0]) +
                         (vertex.y * viewportTransformationMatrix.val[2][1]) +
                         (vertex.z * viewportTransformationMatrix.val[2][2]) +
                         (vertex.t * viewportTransformationMatrix.val[2][3])),
                        vertex.colorId,
                };
                meshVertex.push_back(mshVertex);
                meshVertexWithVp.push_back(mshVertexWithVp);
            }
        }
        allNewVertex.push_back(meshVertex);
        allNewVertexWithVp.push_back(meshVertexWithVp);

    }
    for (int i = 0; i < meshes.size(); ++i) {
        int mesh_type = meshes[i]->type;
        for (int j = 0; j < meshes[i]->numberOfTriangles * 3; j += 3) {
            if (this->cullingEnabled)
                if (backfaceCulling(i, j, camera, allNewVertex) == BACK)
                    continue;

            if (mesh_type == 1) {
                triangleRasterization(i, j, allNewVertexWithVp);
            } else {
                lineRasterization(i, j, meshes[i]->meshId, camera, allNewVertexWithVp);
            }
        }
    }
}

int smallest(int x, int y, int z) {
    int minXY = min(min(x, y), z);
    return minXY < 0 ? 0 : minXY;
}

int largest(int x, int y, int z) {
    int maxYX = max(max(x, y), z);
    return maxYX < 0 ? 0 : maxYX;

}

double f_01(int x, int y, int x_0, int y_0, int x_1, int y_1) {
    return x * (y_0 - y_1) + y * (x_1 - x_0) + x_0 * y_1 - y_0 * x_1;
}

double f_12(int x, int y, int x_1, int y_1, int x_2, int y_2) {
    return x * (y_1 - y_2) + y * (x_2 - x_1) + x_1 * y_2 - y_1 * x_2;
}

double f_20(int x, int y, int x_0, int y_0, int x_2, int y_2) {
    return x * (y_2 - y_0) + y * (x_0 - x_2) + x_2 * y_0 - y_2 * x_0;
}

void Scene::interpolate(int x, int y, const Vec3 &a, const Vec3 &b, const Color &color_a, const Color &color_b) {
    double alpha = 0;
    double dx = abs(b.x - a.x), dy = abs(b.y - a.y);

    dx > dy ? alpha = (x - a.x) / (b.x - a.x) : (y - a.y) / (b.y - a.y);


    double red = (1 - alpha) * (color_a.r) + alpha * color_b.r;
    double green = (1 - alpha) * (color_a.g) + alpha * color_b.g;
    double blue = (1 - alpha) * (color_a.b) + alpha * color_b.b;

    red = makeBetweenZeroAnd255(red);
    green = makeBetweenZeroAnd255(green);
    blue = makeBetweenZeroAnd255(blue);

    Color c = Color(red, green, blue);

    writeToImage(c, x, y);

}

void Scene::writeToImage(const Color &c, int x, int y) {
    image[x][y].r = c.r;
    image[x][y].g = c.g;
    image[x][y].b = c.b;
}

void Scene::triangleRasterization(int i, int j, const vector<vector<Vec3>> &allNewVertexWithVp) {

    double alpha, beta, gamma = 0;
    int x0 = (int) allNewVertexWithVp[i][j].x;
    int x1 = (int) allNewVertexWithVp[i][j + 1].x;
    int x2 = (int) allNewVertexWithVp[i][j + 2].x;

    int y0 = (int) allNewVertexWithVp[i][j].y;
    int y1 = (int) allNewVertexWithVp[i][j + 1].y;
    int y2 = (int) allNewVertexWithVp[i][j + 2].y;

    int x_min = smallest(x0, x1, x2);
    int y_min = smallest(y0, y1, y2);
    int x_max = largest(x0, x1, x2);
    int y_max = largest(y0, y1, y2);

    Color *c0 = this->colorsOfVertices[allNewVertexWithVp[i][j].colorId - 1];
    Color *c1 = this->colorsOfVertices[allNewVertexWithVp[i][j + 1].colorId - 1];
    Color *c2 = this->colorsOfVertices[allNewVertexWithVp[i][j + 2].colorId - 1];

    for (int y = y_min; y <= y_max; y++) {
        for (int x = x_min; x <= x_max; x++) {
            //TODO: Hata görürsen ilk buraya bak!!
            alpha = f_12(x, y, x1, y1, x2, y2) / f_12(x0, y0, x1, y1, x2, y2);
            beta = f_20(x, y, x0, y0, x2, y2) / f_20(x1, y1, x0, y0, x2, y2);
            gamma = f_01(x, y, x0, y0, x1, y1) / f_01(x2, y2, x0, y0, x1, y1);

            if (alpha >= 0 && beta >= 0 && gamma >= 0) {

                auto *color = new Color(alpha * c0->r + beta * c1->r + gamma * c2->r,
                                        alpha * c0->g + beta * c1->g + gamma * c2->g,
                                        alpha * c0->b + beta * c1->b + gamma * c2->b);
                writeToImage(*color, x, y);
            }
        }
    }
}

void Scene::drawLine(const Vec3 &smallerVertex, const Vec3 &biggerVertex, REGION region, const Color &color_a,
                     const Color &color_b) {
    int s_x = (int) smallerVertex.x;
    int b_x = (int) biggerVertex.x;
    int s_y = (int) smallerVertex.y;
    int b_y = (int) biggerVertex.y;

    int x = s_x;
    int y = s_y;
    double t;
    switch (region) {
        case FIRST:
            t = (s_x - b_x) + 0.5 * (s_y - b_y);

            for (y = s_y; y > b_y; --y) {
                interpolate(x, y, smallerVertex, biggerVertex, color_a, color_b);
                if (t < 0) {
                    x += 1;
                    t += (s_x - b_x) + (s_y - b_y);
                } else {
                    t += (s_x - b_x);
                }
            }
            break;
        case SECOND:
            t = (b_y - s_y) + 0.5 * (b_x - s_x);

            for (x = s_x; x < b_x; ++x) {
                interpolate(x, y, smallerVertex, biggerVertex, color_a, color_b);

                if (t < 0) {
                    y -= 1;
                    t += (b_y - s_y) + (b_x - s_x);
                } else {
                    t += (b_y - s_y);
                }
            }
            break;
        case THIRD:
            t = (s_y - b_y) + 0.5 * (b_x - s_x);

            for (x = s_x; x < b_x; ++x) {
                interpolate(x, y, smallerVertex, biggerVertex, color_a, color_b);
                if (t < 0) {
                    y += 1;
                    t += (s_y - b_y) + (b_x - s_x);
                } else {
                    t += (s_y - b_y);
                }
            }
            break;
        case FOURTH:
            t = (s_x - b_x) + 0.5 * (b_y - s_y);

            for (y = s_y; y < b_y; ++y) {
                interpolate(x, y, smallerVertex, biggerVertex, color_a, color_b);

                if (t < 0) {
                    x += 1;
                    t += (s_x - b_x) + (b_y - s_y);
                } else {
                    t += (s_x - b_x);
                }
            }
    }
}

double calculateSlope(const Vec3 &a, const Vec3 &b) {
    return (double) (b.y - a.y) / (double) (b.x - a.x);
}

void Scene::lineRasterization(int i, int j, int id, Camera *cam, vector<vector<Vec3>> vpvertices) {
    Vec3 v0 = vpvertices[i][j];
    Vec3 v1 = vpvertices[i][j + 1];
    Vec3 v2 = vpvertices[i][j + 2];
    double m;

    vector<Vec3> temp = {v0, v1, v2};

    for (int k = 0; k < 3; ++k) {
        Vec3 a = temp[k];
        Vec3 b = temp[(k + 1) % 3];
        int x, y, d;

        Vec3 aclipped, bclipped;
        Color color_a(*(this->colorsOfVertices[a.colorId - 1])), color_b(*(this->colorsOfVertices[b.colorId - 1]));
        bool isVisible = LiangBarskyAlgorithm(a, b, cam, aclipped, bclipped, color_a, color_b, *this);
        if (!isVisible) continue;
        a = aclipped;
        b = bclipped;


        //mandatory to get rid of clipping errors!!
        if ((int) a.x > cam->verRes - 1) a.x = cam->verRes - 1;
        if ((int) a.y > cam->horRes - 1) a.x = cam->horRes - 1;
        if ((int) a.x < 0) a.x = 0;
        if ((int) a.y < 0) a.y = 0;

        if ((int) b.x > cam->verRes - 1) b.x = cam->verRes - 1;
        if ((int) b.y > cam->horRes - 1) b.x = cam->horRes - 1;
        if ((int) b.x < 0) b.x = 0;
        if ((int) b.y < 0) b.y = 0;


        m = calculateSlope(a, b);

        if (m < -1)
            a.x > b.x ? drawLine(b, a, FIRST, color_b, color_a) : drawLine(a, b, FIRST, color_a, color_b);
        else if (m < 0)
            a.x > b.x ? drawLine(b, a, SECOND, color_b, color_a) : drawLine(a, b, SECOND, color_a, color_b);
        else if (m < 1)
            a.x > b.x ? drawLine(b, a, THIRD, color_b, color_a) : drawLine(a, b, THIRD, color_a, color_b);
        else
            a.x > b.x ? drawLine(b, a, FOURTH, color_b, color_a) : drawLine(a, b, FOURTH, color_a, color_b);


    }
}


Matrix4 Scene::ModelingTransformation(Mesh *mesh) {
    Matrix4 modeledMatrix = getIdentityMatrix();
    Matrix4 preComputedMatrix = getIdentityMatrix();

    for (int i = 0; i < mesh->numberOfTransformations; i++) {
        switch (mesh->transformationTypes[i]) {
            case 'r':
                preComputedMatrix = computeRotationMatrix(this->rotations[mesh->transformationIds[i] - 1]);
                modeledMatrix = multiplyMatrixWithMatrix(preComputedMatrix, modeledMatrix);
                break;
            case 't':
                preComputedMatrix = computeTranslationMatrix(
                        this->translations[mesh->transformationIds[i] - 1]);
                modeledMatrix = multiplyMatrixWithMatrix(preComputedMatrix, modeledMatrix);
                break;
            case 's':
                preComputedMatrix = computeScalingMatrix(this->scalings[mesh->transformationIds[i] - 1]);
                modeledMatrix = multiplyMatrixWithMatrix(preComputedMatrix, modeledMatrix);
                break;
        }
    }
    return modeledMatrix;
}

/*
    Parses XML file
*/
Scene::Scene(
        const char *xmlPath) {
    const char *str;
    XMLDocument xmlDoc;
    XMLElement *pElement;

    xmlDoc.LoadFile(xmlPath);

    XMLNode *pRoot = xmlDoc.FirstChild();

    // read background color
    pElement = pRoot->FirstChildElement("BackgroundColor");
    str = pElement->GetText();
    sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

    // read backfaceCulling
    pElement = pRoot->FirstChildElement("Culling");
    if (pElement != NULL) {
        str = pElement->GetText();

        if (strcmp(str, "enabled") == 0) {
            cullingEnabled = true;
        } else {
            cullingEnabled = false;
        }
    }

    // read cameras
    pElement = pRoot->FirstChildElement("Cameras");
    XMLElement *pCamera = pElement->FirstChildElement("Camera");
    XMLElement *camElement;
    while (pCamera != NULL) {
        Camera *cam = new Camera();

        pCamera->QueryIntAttribute("id", &cam->cameraId);

        // read projection type
        str = pCamera->Attribute("type");

        if (strcmp(str, "orthographic") == 0) {
            cam->projectionType = 0;
        } else {
            cam->projectionType = 1;
        }

        camElement = pCamera->FirstChildElement("Position");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

        camElement = pCamera->FirstChildElement("Gaze");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

        camElement = pCamera->FirstChildElement("Up");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

        cam->gaze = normalizeVec3(cam->gaze);
        cam->u = crossProductVec3(cam->gaze, cam->v);
        cam->u = normalizeVec3(cam->u);

        cam->w = inverseVec3(cam->gaze);
        cam->v = crossProductVec3(cam->u, cam->gaze);
        cam->v = normalizeVec3(cam->v);

        camElement = pCamera->FirstChildElement("ImagePlane");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
               &cam->left, &cam->right, &cam->bottom, &cam->top,
               &cam->near, &cam->far, &cam->horRes, &cam->verRes);

        camElement = pCamera->FirstChildElement("OutputName");
        str = camElement->GetText();
        cam->outputFileName = string(str);

        cameras.push_back(cam);

        pCamera = pCamera->NextSiblingElement("Camera");
    }

    // read vertices
    pElement = pRoot->FirstChildElement("Vertices");
    XMLElement *pVertex = pElement->FirstChildElement("Vertex");
    int vertexId = 1;

    while (pVertex != NULL) {
        Vec3 *vertex = new Vec3();
        Color *color = new Color();

        vertex->colorId = vertexId;

        str = pVertex->Attribute("position");
        sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

        str = pVertex->Attribute("color");
        sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

        vertices.push_back(vertex);
        colorsOfVertices.push_back(color);

        pVertex = pVertex->NextSiblingElement("Vertex");

        vertexId++;
    }

    // read translations
    pElement = pRoot->FirstChildElement("Translations");
    XMLElement *pTranslation = pElement->FirstChildElement("Translation");
    while (pTranslation != NULL) {
        Translation *translation = new Translation();

        pTranslation->QueryIntAttribute("id", &translation->translationId);

        str = pTranslation->Attribute("value");
        sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

        translations.push_back(translation);

        pTranslation = pTranslation->NextSiblingElement("Translation");
    }

    // read scalings
    pElement = pRoot->FirstChildElement("Scalings");
    XMLElement *pScaling = pElement->FirstChildElement("Scaling");
    while (pScaling != NULL) {
        Scaling *scaling = new Scaling();

        pScaling->QueryIntAttribute("id", &scaling->scalingId);
        str = pScaling->Attribute("value");
        sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

        scalings.push_back(scaling);

        pScaling = pScaling->NextSiblingElement("Scaling");
    }

    // read rotations
    pElement = pRoot->FirstChildElement("Rotations");
    XMLElement *pRotation = pElement->FirstChildElement("Rotation");
    while (pRotation != NULL) {
        Rotation *rotation = new Rotation();

        pRotation->QueryIntAttribute("id", &rotation->rotationId);
        str = pRotation->Attribute("value");
        sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

        rotations.push_back(rotation);

        pRotation = pRotation->NextSiblingElement("Rotation");
    }

    // read meshes
    pElement = pRoot->FirstChildElement("Meshes");

    XMLElement *pMesh = pElement->FirstChildElement("Mesh");
    XMLElement *meshElement;
    while (pMesh != NULL) {
        Mesh *mesh = new Mesh();

        pMesh->QueryIntAttribute("id", &mesh->meshId);

        // read projection type
        str = pMesh->Attribute("type");

        if (strcmp(str, "wireframe") == 0) {
            mesh->type = 0;
        } else {
            mesh->type = 1;
        }

        // read mesh transformations
        XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
        XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

        while (pTransformation != NULL) {
            char transformationType;
            int transformationId;

            str = pTransformation->GetText();
            sscanf(str, "%c %d", &transformationType, &transformationId);

            mesh->transformationTypes.push_back(transformationType);
            mesh->transformationIds.push_back(transformationId);

            pTransformation = pTransformation->NextSiblingElement("Transformation");
        }

        mesh->numberOfTransformations = mesh->transformationIds.size();

        // read mesh faces
        char *row;
        char *clone_str;
        int v1, v2, v3;
        XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
        clone_str = strdup(str);

        row = strtok(clone_str, "\n");
        while (row != NULL) {
            int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

            if (result != EOF) {
                mesh->triangles.push_back(Triangle(v1, v2, v3));
            }
            row = strtok(NULL, "\n");
        }
        mesh->numberOfTriangles = mesh->triangles.size();
        meshes.push_back(mesh);

        pMesh = pMesh->NextSiblingElement("Mesh");
    }
}

/*
    Initializes image with background color
*/
void Scene::initializeImage(Camera *camera) {
    if (this->image.empty()) {
        for (int i = 0; i < camera->horRes; i++) {
            vector<Color> rowOfColors;

            for (int j = 0; j < camera->verRes; j++) {
                rowOfColors.push_back(this->backgroundColor);
            }

            this->image.push_back(rowOfColors);
        }
    } else {
        for (int i = 0; i < camera->horRes; i++) {
            for (int j = 0; j < camera->verRes; j++) {
                this->image[i][j].r = this->backgroundColor.r;
                this->image[i][j].g = this->backgroundColor.g;
                this->image[i][j].b = this->backgroundColor.b;
            }
        }
    }
}

/*
    If given value is less than 0, converts value to 0.
    If given value is more than 255, converts value to 255.
    Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value) {
    if (value >= 255.0)
        return 255;
    if (value <= 0.0)
        return 0;
    return (int) (value);
}

/*
    Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera) {
    ofstream fout;

    fout.open(camera->outputFileName.c_str());

    fout << "P3" << endl;
    fout << "# " << camera->outputFileName << endl;
    fout << camera->horRes << " " << camera->verRes << endl;
    fout << "255" << endl;

    for (int j = camera->verRes - 1; j >= 0; j--) {
        for (int i = 0; i < camera->horRes; i++) {
            fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
        }
        fout << endl;
    }
    fout.close();
}

/*
    Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
    os_type == 1 		-> Ubuntu
    os_type == 2 		-> Windows
    os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType) {
    string command;

    // call command on Ubuntu
    if (osType == 1) {
        command = "convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

        // call command on Windows
    else if (osType == 2) {
        command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    } else if (osType == 3) {
        command = "pnmtopng " + ppmFileName + " > " + ppmFileName + ".png";
        system(command.c_str());
    }
        // default action - don't do conversion
    else {

    }
}