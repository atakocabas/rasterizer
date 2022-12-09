#include <iostream>
#include <iomanip>
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

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/
Matrix4 ModelingTransformation(Mesh *mesh);


void Scene::forwardRenderingPipeline(Camera *camera) {
    // TODO: Implement this function.
    /*  Matrix4 viewportMatrix = computeViewingTransformationMatrix(camera);

      if (camera->projectionType == 0) {//ORTHOGRAPHIC
          Matrix4 orthographicProjectionMatrix = calculateOrthographicProjection(camera);
          camViewMatrix = multiplyMatrixWithMatrix(camViewMatrix, orthographicProjectionMatrix);
      } else {
          Matrix4 perspectiveProjectionMatrix = calculatePerspectiveProjection(camera);
          camViewMatrix = multiplyMatrixWithMatrix(camViewMatrix, perspectiveProjectionMatrix);
      }
      camViewMatrix = multiplyMatrixWithMatrix(viewportMatrix, camViewMatrix);
      */


    vector<vector<Vec3>> allNewVertex;
    vector<Vec3> meshVertex;
    Matrix4 orthPerspectiveProjectionMatrix;
    Matrix4 camViewMatrix = computeViewingTransformationMatrix(camera);
    if (camera->projectionType == 0) {//ORTHOGRAPHIC
         orthPerspectiveProjectionMatrix = calculateOrthographicProjection(camera);
        orthPerspectiveProjectionMatrix = multiplyMatrixWithMatrix(camViewMatrix, orthPerspectiveProjectionMatrix);
    } else {
        Matrix4 perspectiveProjectionMatrix = calculatePerspectiveProjection(camera);
        orthPerspectiveProjectionMatrix = multiplyMatrixWithMatrix(camViewMatrix, perspectiveProjectionMatrix);
    }
    Matrix4 viewportTransformationMatrix = calculateViewportMatrix(camera);
    for (auto &mesh: this->meshes) {
        Matrix4 modelingViewMatrix = ModelingTransformation(mesh);
        modelingViewMatrix = multiplyMatrixWithMatrix(orthPerspectiveProjectionMatrix, multiplyMatrixWithMatrix(camViewMatrix,modelingViewMatrix));

        for (auto &triangle: mesh->triangles) {

            //TODO: Check Ids if any range error thrown
            Vec3 *new_vertex = this->vertices[triangle.getFirstVertexId() - 1];

            Vec4 vertex = {
                    new_vertex->x, new_vertex->y, new_vertex->z, 1, new_vertex->colorId,
            };
            vertex = multiplyMatrixWithVec4(modelingViewMatrix, vertex);
            Vec3 mshVertex = {
                    new_vertex->x, new_vertex->y, new_vertex->z, new_vertex->colorId,
            };
            meshVertex.push_back(mshVertex);
        }
        allNewVertex.push_back(meshVertex);

    }
}


Matrix4 Scene::ModelingTransformation(Mesh *mesh) {
    Matrix4 modeledMatrix = getIdentityMatrix();
    Matrix4 preComputedMatrix = getIdentityMatrix();

    for (int i = 0; i < mesh->numberOfTransformations; i++) {
        switch (mesh->transformationTypes[i]) {
            case 'r':
                preComputedMatrix = computeRotationMatrix(this->rotations[mesh->transformationIds[i] - 1]);
                modeledMatrix = multiplyMatrixWithMatrix(preComputedMatrix,modeledMatrix );
                break;
            case 't':
                preComputedMatrix = computeTranslationMatrix(
                        this->translations[mesh->transformationIds[i] - 1]);
                modeledMatrix = multiplyMatrixWithMatrix(preComputedMatrix,modeledMatrix );
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

    // read culling
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
    }

        // default action - don't do conversion
    else {
    }
}