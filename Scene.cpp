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

void Scene::forwardRenderingPipeline(Camera *camera)
{
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
    vector<vector<Vec3>> allNewVertexWithVp;
    vector<Vec3> meshVertex;
    vector<Vec3> meshVertexWithVp;
    Matrix4 orthPerspectiveProjectionMatrix;
    Matrix4 camViewMatrix = computeViewingTransformationMatrix(camera);
    if (camera->projectionType == 0)
    { // ORTHOGRAPHIC
        orthPerspectiveProjectionMatrix = calculateOrthographicProjection(camera);
    }
    else
    {
        orthPerspectiveProjectionMatrix = calculatePerspectiveProjection(camera);
    }
    Matrix4 viewportTransformationMatrix = calculateViewportMatrix(camera);
    for (auto &mesh : this->meshes)
    {
        Matrix4 modelingViewMatrix = ModelingTransformation(mesh);
        modelingViewMatrix = multiplyMatrixWithMatrix(orthPerspectiveProjectionMatrix,
                                                      multiplyMatrixWithMatrix(camViewMatrix, modelingViewMatrix));

        for (auto &triangle : mesh->triangles)
        {
            for (auto &vertex_id : triangle.vertexIds)
            {
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

                int k, l;
                double total;

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
                // ViewPort here?
                meshVertex.push_back(mshVertex);
                meshVertexWithVp.push_back(mshVertexWithVp);
            }
        }
        allNewVertex.push_back(meshVertex);
        allNewVertexWithVp.push_back(meshVertexWithVp);

        // ROUNDING to set its pixel in vp
    }
    bool checkCulling = this->cullingEnabled;
    for (int i = 0; i < meshes.size(); ++i)
    {
        int mesh_type = meshes[i]->type;
        for (int j = 0; j < meshes[i]->numberOfTriangles * 3; j += 3)
        {
            if (checkCulling && culling(i, j, camera, allNewVertex))
                continue;

            // MID POINT            // RASTERIZATION
            if (mesh_type == 1)
            {
                rasterization(i, j, allNewVertexWithVp);
            }
            else
            {
                midPoint(i, j, meshes[i]->meshId, camera, allNewVertexWithVp);
            }
        }
    }
}

int smallest(int x, int y, int z)
{
    return min(min(x, y), z);
}

int largest(int x, int y, int z)
{
    return max(max(x, y), z);
}

int f_01(int x, int y, int x_0, int y_0, int x_1, int y_1)
{
    return x * (y_0 - y_1) + y * (x_1 - x_0) + x_0 * y_1 - y_0 * x_1;
}

int f_12(int x, int y, int x_1, int y_1, int x_2, int y_2)
{
    return x * (y_1 - y_2) + y * (x_2 - x_1) + x_1 * y_2 - y_1 * x_2;
}

int f_20(int x, int y, int x_0, int y_0, int x_2, int y_2)
{
    return x * (y_2 - y_0) + y * (x_0 - x_2) + x_2 * y_0 - y_2 * x_0;
}
/*
 * void drawLine(Point p1, Point p2, Color color) {
  int dx = p2.x - p1.x;
  int dy = p2.y - p1.y;

  // Calculate the slope of the line
  float m = (float)dy / dx;

  // Check if the line is steep
  if (std::abs(m) > 1) {
    // If the line is steep, swap the x and y coordinates
    std::swap(p1.x, p1.y);
    std::swap(p2.x, p2.y);
    std::swap(dx, dy);

    // Recalculate the slope
    m = (float)dy / dx;
  }

  // Check if the line has a negative slope
  if (p1.x > p2.x) {
    // If the line has a negative slope, swap the points
    std::swap(p1, p2);
  }

  // Calculate the midpoint of the line
  int d = dy - m * dx;
  int y = p1.y;

  // Loop over the x-coordinates of the line
  for (int x = p1.x; x <= p2.x; x++) {
    // Check if the line is steep
    if (std::abs(m) > 1) {
      // If the line is steep, swap the x and y coordinates
      setPixelColor(y, x, color);
    } else {
      // If the line is not steep, draw the pixel with the original coordinates
      setPixelColor(x, y, color);
    }

    // Check if the midpoint is above or below the line
    if (d < 0) {
      // If the midpoint is below the line, move to the next pixel in the y-direction
      y++;
      d += dy;
    } else {
      // If the midpoint is above or on the line, move to the next pixel in the x-direction
      d -= dx;
    }
  }
}
 */
void Scene::rasterization(int i, int j, vector<vector<Vec3>> allNewVertexWithVp)
{

    int x_0 = allNewVertexWithVp[i][j].x;
    int x_1 = allNewVertexWithVp[i][j + 1].x;
    int x_2 = allNewVertexWithVp[i][j + 2].x;

    int y_0 = allNewVertexWithVp[i][j].y;
    int y_1 = allNewVertexWithVp[i][j + 1].y;
    int y_2 = allNewVertexWithVp[i][j + 2].y;

    Color *c0 = this->colorsOfVertices[allNewVertexWithVp[i][j].colorId - 1];
    Color *c1 = this->colorsOfVertices[allNewVertexWithVp[i][j + 1].colorId - 1];
    Color *c2 = this->colorsOfVertices[allNewVertexWithVp[i][j + 2].colorId - 1];

    int x_min = smallest(x_0, x_1, x_2);
    int y_min = smallest(y_0, y_1, y_2);

    int x_max = largest(x_0, x_1, x_2);
    int y_max = largest(y_0, y_1, y_2);

    for (int y = y_min; y <= y_max; y++)
    {
        for (int x = x_min; x <= x_max; x++)
        {
            double alpha = (double)f_12(x, y, x_1, y_1, x_2, y_2) / (double)f_12(x_0, y_0, x_1, y_1, x_2, y_2);
            double beta = (double)f_20(x, y, x_0, y_0, x_2, y_2) / (double)f_20(x_1, y_1, x_0, y_0, x_2, y_2);
            double gama = (double)f_01(x, y, x_0, y_0, x_1, y_1) / (double)f_01(x_2, y_2, x_0, y_0, x_1, y_1);

            if (alpha >= 0 && beta >= 0 && gama >= 0)
            {
                image[x][y].r = alpha * c0->r + beta * c1->r + gama * c2->r;
                image[x][y].g = alpha * c0->g + beta * c1->g + gama * c2->g;
                image[x][y].b = alpha * c0->b + beta * c1->b + gama * c2->b;
            }
        }
    }
}

void Scene::midPoint(int i, int j, int id, Camera *cam, vector<vector<Vec3>> vpvertices)
{
    Vec3 v0 = vpvertices[i][j];
    Vec3 v1 = vpvertices[i][j + 1];
    Vec3 v2 = vpvertices[i][j + 2];
    double m;
    Camera c = *cam;

    vector<Vec3> temp = {v0, v1, v2};

    for (int k = 0; k < 3; ++k)
    {
        Vec3 a = temp[k];
        Vec3 b = temp[(k + 1) % 3];

        Vec3 aclipped, bclipped;
        LiangBarskyAlgorithm(a, b, cam, aclipped, bclipped);
        if (areEqualVec3(aclipped, bclipped))
            continue;
        a = aclipped;
        b = bclipped;
        m = calculateSlope(a, b);
        if (m == 0)
            continue;
        if ((b.y - a.y) >= 0 && (b.x - a.x) >= 0)
		{ // 1.çeyrek
			if (m > 0 && m <= 1)
			{ //1.bölge
				int x = a.x;
				int y = a.y;
				// y=a.y;
				double M = 1.0 * (a.y - b.y) + 0.5 * (b.x - a.x);
				for (x = a.x; x < int(b.x) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; x = x + 1)
				{
					draw(x, y, a, b);
					if (M < 0)
					{ // NE
						y += 1;
						M += ((a.y - b.y) + (b.x - a.x));
					}
					else
					{ // E
						M += (a.y - b.y);
					}
				}
			}
			else if (m > 1)
			{ //2.bölge
				int x = a.x;
				int y = a.y;
				double M = 0.5 * (a.y - b.y) + 1.0 * (b.x - a.x);
				for (y = a.y; y < int(b.y) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; y = y + 1)
				{
					draw(x, y, a, b);
					if (M <= 0)
					{ // N
						M += (b.x - a.x);
					}
					else
					{ // NE
						x += 1;
						M += ((a.y - b.y) + (b.x - a.x));
					}
				}
			}
		}
		else if ((b.y - a.y) > 0 && (b.x - a.x) < 0) // 2.çeyrek
		{
			if (m < -1)
			{ //2.bölge
				int x = a.x;
				int y = a.y;
				double M = -0.5 * (a.y - b.y) + 1.0 * (b.x - a.x);
				for (y = a.y; y < int(b.y) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; y = y + 1)
				{
					draw(x, y, a, b);
					if (M > 0)
					{ // N
						M += ((b.x - a.x));
					}
					else
					{ // NW
						x -= 1;
						M += (-(a.y - b.y) + (b.x - a.x));
					}
				}
			}
			else if (m < 0 && m >= -1)
			{ // 1.bölge
				int x = a.x;
				int y = a.y;
				double M = -1.0 * (a.y - b.y) + 0.5 * (b.x - a.x);
				for (x = a.x; x > int(b.x) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; x = x - 1)
				{
					draw(x, y, a, b);
					if (M > 0)
					{ // NW
						y += 1;
						M += (-(a.y - b.y) + (b.x - a.x));
					}
					else
					{ // W
						M += (-(a.y - b.y));
					}
				}
			}
		}
		else if ((b.y - a.y) <= 0 && (b.x - a.x) <= 0) // 3.çeyrek
		{
			if (m > 1)
			{ // 2.bölge
				int x = a.x;
				int y = a.y;
				double M = -0.5 * (a.y - b.y) - 1.0 * (b.x - a.x);
				for (y = a.y; y > int(b.y) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; y = y - 1)
				{
					draw(x, y, a, b);
					if (M > 0)
					{ // SW
						x -= 1;
						M += (-(a.y - b.y) - (b.x - a.x));
					}
					else
					{ // S
						M += (-(b.x - a.x));
					}
				}
			}
			else if (m > 0 && m <= 1)
			{ //1.bölge
				int x = a.x;
				int y = a.y;
				double M = -1.0 * (a.y - b.y) - 0.5 * (b.x - a.x);
				for (x = a.x; x > int(b.x) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; x = x - 1)
				{
					draw(x, y, a, b);
					if (M > 0)
					{ // W
						M += (-(a.y - b.y));
					}
					else
					{ // SW
						y -= 1;
						M += (-(a.y - b.y) - (b.x - a.x));
					}
				}
			}
		}
		else if ((b.y - a.y) < 0 && (b.x - a.x) > 0) //4.çeyrek
		{
			if (m < 0 && m >= -1)
			{ //1.bölge
				int x = a.x;
				int y = a.y;
				double M = 1.0 * (a.y - b.y) - 0.5 * (b.x - a.x);
				for (x = a.x; x < int(b.x) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; x = x + 1)
				{
					draw(x, y, a, b);
					if (M <= 0)
					{ // E
						M += ((a.y - b.y));
					}
					else
					{ // SE
						y -= 1;
						M += ((a.y - b.y) - (b.x - a.x));
					}
				}
			}
			else if (m < -1)
			{ //2.bölge
				int x = a.x;
				int y = a.y;
				double M = 0.5 * (a.y - b.y) - 1.0 * (b.x - a.x);
				for (y = a.y; y > int(b.y) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; y = y - 1)
				{
					draw(x, y, a, b);
					if (M < 0)
					{ // SE
						x += 1;
						M += ((a.y - b.y) - (b.x - a.x));
					}
					else
					{ // S
						M += (-(b.x - a.x));
					}
				}
			}
		}
    }
}

void Scene::draw(int x, int y, Vec3 a, Vec3 b)
{
    double alphaX = (x - a.x) / (b.x - a.x);
    double alphaY = (y - a.y) / (b.y - a.y);

    Color *color_a = colorsOfVertices[a.colorId - 1];
    Color *color_b = colorsOfVertices[b.colorId - 1];
    double cX_r = (1 - alphaX) * (color_a->r) + alphaX * color_b->r;
    double cX_g = (1 - alphaX) * (color_a->g) + alphaX * color_b->g;
    double cX_b = (1 - alphaX) * (color_a->b) + alphaX * color_b->b;

    if (cX_r > 255)
        cX_r = 255;
    else if (cX_r < 0)
        cX_r = 0;
    if (cX_g > 255)
        cX_g = 255;
    else if (cX_g < 0)
        cX_g = 0;
    if (cX_b > 255)
        cX_b = 255;
    else if (cX_b < 0)
        cX_b = 0;

    Color c = Color(cX_r, cX_g, cX_b);

    image[x][y].r = c.r;
    image[x][y].g = c.g;
    image[x][y].b = c.b;
}

Matrix4 Scene::ModelingTransformation(Mesh *mesh)
{
    Matrix4 modeledMatrix = getIdentityMatrix();
    Matrix4 preComputedMatrix = getIdentityMatrix();

    for (int i = 0; i < mesh->numberOfTransformations; i++)
    {
        switch (mesh->transformationTypes[i])
        {
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
    const char *xmlPath)
{
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
    if (pElement != NULL)
    {
        str = pElement->GetText();

        if (strcmp(str, "enabled") == 0)
        {
            cullingEnabled = true;
        }
        else
        {
            cullingEnabled = false;
        }
    }

    // read cameras
    pElement = pRoot->FirstChildElement("Cameras");
    XMLElement *pCamera = pElement->FirstChildElement("Camera");
    XMLElement *camElement;
    while (pCamera != NULL)
    {
        Camera *cam = new Camera();

        pCamera->QueryIntAttribute("id", &cam->cameraId);

        // read projection type
        str = pCamera->Attribute("type");

        if (strcmp(str, "orthographic") == 0)
        {
            cam->projectionType = 0;
        }
        else
        {
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

    while (pVertex != NULL)
    {
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
    while (pTranslation != NULL)
    {
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
    while (pScaling != NULL)
    {
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
    while (pRotation != NULL)
    {
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
    while (pMesh != NULL)
    {
        Mesh *mesh = new Mesh();

        pMesh->QueryIntAttribute("id", &mesh->meshId);

        // read projection type
        str = pMesh->Attribute("type");

        if (strcmp(str, "wireframe") == 0)
        {
            mesh->type = 0;
        }
        else
        {
            mesh->type = 1;
        }

        // read mesh transformations
        XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
        XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

        while (pTransformation != NULL)
        {
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
        while (row != NULL)
        {
            int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

            if (result != EOF)
            {
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
void Scene::initializeImage(Camera *camera)
{
    if (this->image.empty())
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            vector<Color> rowOfColors;

            for (int j = 0; j < camera->verRes; j++)
            {
                rowOfColors.push_back(this->backgroundColor);
            }

            this->image.push_back(rowOfColors);
        }
    }
    else
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            for (int j = 0; j < camera->verRes; j++)
            {
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
int Scene::makeBetweenZeroAnd255(double value)
{
    if (value >= 255.0)
        return 255;
    if (value <= 0.0)
        return 0;
    return (int)(value);
}

/*
    Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
    ofstream fout;

    fout.open(camera->outputFileName.c_str());

    fout << "P3" << endl;
    fout << "# " << camera->outputFileName << endl;
    fout << camera->horRes << " " << camera->verRes << endl;
    fout << "255" << endl;

    for (int j = camera->verRes - 1; j >= 0; j--)
    {
        for (int i = 0; i < camera->horRes; i++)
        {
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
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
    string command;

    // call command on Ubuntu
    if (osType == 1)
    {
        command = "convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

    // call command on Windows
    else if (osType == 2)
    {
        command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

    // default action - don't do conversion
    else
    {
    }
}