#include <iostream>
#include <cmath>
#include "Helpers.h"
#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Translation.h"
#include "Scaling.h"
#include "Rotation.h"
#include <vector>

using namespace std;

ORDER backfaceCulling(int i, int j, Camera *cam, vector<vector<Vec3>> cameraview) {

    Vec3 v0 = cameraview[i][j];
    Vec3 v1 = cameraview[i][j + 1];
    Vec3 v2 = cameraview[i][j + 2];

    Vec3 eye = subtractVec3(cam->pos, v0);
    Vec3 normal = crossProductVec3(subtractVec3(v1, v0), subtractVec3(v2, v0));

    if (dotProductVec3(normal, v1) > 0)
        return FRONT;
    else
        return BACK;
}

Matrix4 calculateViewportMatrix(Camera *camera) {
    Matrix4 res = getIdentityMatrix();

    res.val[0][0] = camera->horRes / 2.0;
    res.val[0][3] = ((camera->horRes - 1) / 2.0);
    res.val[1][1] = camera->verRes / 2.0;
    res.val[1][3] = ((camera->verRes - 1) / 2.0);
    res.val[2][2] = 0.5;
    res.val[2][3] = 0.5;

    res.val[3][3] = 0;

    return res;
}

bool LiangBarskyAlgorithm(const Vec3 &v0, const Vec3 &v1, Camera *cam, Vec3 &v0new, Vec3 &v1new, Color &color_a,
                          Color &color_b,
                          const Scene &scene) {
    Color v0color = *(scene.colorsOfVertices.at(v0.colorId - 1)), v1color = *(scene.colorsOfVertices.at(
            v1.colorId - 1));
    double te = 0.0, tl = 1.0;
    double dx = v1.x - v0.x, dy = v1.y - v0.y;
    bool visible = false;
    double alpha;
    v0new.x = v0.x;
    v0new.y = v0.y;
    v0new.colorId = v0.colorId;

    v1new.x = v1.x;
    v1new.y = v1.y;
    v1new.colorId = v1.colorId;


    if (isVisible(dx, 0 - v0.x, te, tl))
        if (isVisible(-dx, v0.x - cam->verRes, te, tl))
            if (isVisible(dy, 0 - v0.y, te, tl))
                if (isVisible(-dy, v0.y - cam->horRes, te, tl)) {
                    visible = true;
                    if (tl < 1) {
                        v1new.x = v0.x + dx * tl;
                        v1new.y = v0.y + dy * tl;
                        alpha = (v1new.x - v0.x) / (v1.x - v0.x);
                        color_b.r = (1 - alpha) * v0color.r + alpha * v1color.r;
                        color_b.g = (1 - alpha) * v0color.g + alpha * v1color.g;
                        color_b.b = (1 - alpha) * v0color.b + alpha * v1color.b;
                    }
                    if (te > 0) {
                        v0new.x = v0.x + dx * te;
                        v0new.y = v0.y + dy * te;
                        alpha = (v0new.x - v0.x) / (v1.x - v0.x);
                        color_a.r = (1 - alpha) * v0color.r + alpha * v1color.r;
                        color_a.g = (1 - alpha) * v0color.g + alpha * v1color.g;
                        color_a.b = (1 - alpha) * v0color.b + alpha * v1color.b;
                    }
                }
    return visible;
}

bool isVisible(double den, double num, double &te, double &tl) {
    double t;
    if (den > 0) {
        t = num / den;
        if (t > tl)
            return false;
        if (t > te)
            te = t;
    } else if (den < 0) {
        t = num / den;
        if (t < te)
            return false;
        if (t < tl)
            tl = t;
    } else if (num > 0)
        return false; // parallel line
    return true;
}

Matrix4 calculatePerspectiveProjection(Camera *camera) {
    Matrix4 p2o = getIdentityMatrix();
    Matrix4 ortoMatrix = calculateOrthographicProjection(camera);
    p2o.val[0][0] = camera->near;
    p2o.val[1][1] = camera->near;

    p2o.val[2][2] = camera->far + camera->near;
    p2o.val[2][3] = camera->far * camera->near;

    p2o.val[3][2] = -1;
    p2o.val[3][3] = 0;

    return multiplyMatrixWithMatrix(ortoMatrix, p2o);
}

Matrix4 calculateOrthographicProjection(Camera *camera) {
    Matrix4 orthographic_matrix = getIdentityMatrix();

    // DONE: What if right-left = 0 -> Of course the plane would be a line if it happens!
    orthographic_matrix.val[0][0] = 2 / (camera->right - camera->left);
    orthographic_matrix.val[0][3] = -(camera->right + camera->left) / (camera->right - camera->left);

    orthographic_matrix.val[1][1] = 2 / (camera->top - camera->bottom);
    orthographic_matrix.val[1][3] = -(camera->top + camera->bottom) / (camera->top - camera->bottom);

    orthographic_matrix.val[2][2] = -2 / (camera->far - camera->near);
    orthographic_matrix.val[2][3] = -(camera->far + camera->near) / (camera->far - camera->near);

    return orthographic_matrix;
}

Matrix4 computeViewingTransformationMatrix(Camera *cam) {
    Matrix4 res;
    Matrix4 t = getIdentityMatrix(), m = getIdentityMatrix();
    t.val[0][3] = -cam->pos.x;
    t.val[1][3] = -cam->pos.y;
    t.val[2][3] = -cam->pos.z;

    m.val[0][0] = cam->u.x;
    m.val[0][1] = cam->u.y;
    m.val[0][2] = cam->u.z;
    m.val[1][0] = cam->v.x;
    m.val[1][1] = cam->v.y;
    m.val[1][2] = cam->v.z;
    m.val[2][0] = cam->w.x;
    m.val[2][1] = cam->w.y;
    m.val[2][2] = cam->w.z;

    res = multiplyMatrixWithMatrix(m, t);
    return res;
}

Matrix4 computeTranslationMatrix(Translation *t) {
    Matrix4 m = getIdentityMatrix();
    m.val[0][3] = t->tx;
    m.val[1][3] = t->ty;
    m.val[2][3] = t->tz;

    return m;
}

Matrix4 computeRotationMatrix(Rotation *rotation) {
    Rotation *r = rotation;
    double angle = r->angle;
    double theta = angle * M_PI / 180;

    Vec3 u(r->ux, r->uy, r->uz, 0);
    Matrix4 res;

    double temp_min = abs(u.x);
    Vec3 v = Vec3(0, -u.z, u.y, -1);

    if (abs(u.y) < temp_min) {
        temp_min = abs(u.y);
        v.x = -u.z;
        v.y = 0;
        v.z = u.x;
    } else if (abs(u.z) < temp_min) {
        v.x = -u.y;
        v.y = u.x;
        v.z = 0;
    }
    v = normalizeVec3(v);
    u = normalizeVec3(u);
    Vec3 w = crossProductVec3(u, v);


    double M[4][4] = {{u.x, u.y, u.z, 0},
                      {v.x, v.y, v.z, 0},
                      {w.x, w.y, w.z, 0},
                      {0,   0,   0,   1}};

    double M_INVERSE[4][4] = {{u.x, v.x, w.x, 0},
                              {u.y, v.y, w.y, 0},
                              {u.z, v.z, w.z, 0},
                              {0,   0,   0,   1}};

    double Rotation[4][4] = {{1, 0,          0,           0},
                             {0, cos(theta), -sin(theta), 0},
                             {0, sin(theta), cos(theta),  0},
                             {0, 0,          0,           1}};

    Matrix4 tmp;
    tmp = multiplyMatrixWithMatrix(M_INVERSE, Rotation);
    res = multiplyMatrixWithMatrix(tmp, M);

    return res;
}

Matrix4 computeScalingMatrix(Scaling *s) {
    Matrix4 m = getIdentityMatrix();
    m.val[0][0] = s->sx;
    m.val[1][1] = s->sy;
    m.val[2][2] = s->sz;

    return m;
}

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(const Vec3 &a, const Vec3 &b) {
    Vec3 result;

    result.x = a.y * b.z - b.y * a.z;
    result.y = b.x * a.z - a.x * b.z;
    result.z = a.x * b.y - b.x * a.y;

    return result;
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(const Vec3 &a, const Vec3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(const Vec3 &v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(const Vec3 &v) {
    Vec3 result;
    double d;

    d = magnitudeOfVec3(v);
    result.x = v.x / d;
    result.y = v.y / d;
    result.z = v.z / d;

    return result;
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(const Vec3 &v) {
    Vec3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;

    return result;
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(const Vec3 &a, const Vec3 &b) {
    Vec3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(const Vec3 &a, const Vec3 &b) {
    Vec3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(const Vec3 &v, double c) {
    Vec3 result;
    result.x = v.x * c;
    result.y = v.y * c;
    result.z = v.z * c;

    return result;
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v) {
    cout << "(" << v.x << "," << v.y << "," << v.z << ")" << endl;
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b) {

    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((ABS((a.x - b.x)) < EPSILON) && (ABS((a.y - b.y)) < EPSILON) && (ABS((a.z - b.z)) < EPSILON)) {
        return 1;
    } else {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
 */
Matrix4 getIdentityMatrix() {
    Matrix4 result;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) {
                result.val[i][j] = 1.0;
            } else {
                result.val[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(const Matrix4 &m1, const Matrix4 &m2) {
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            total = 0;
            for (int k = 0; k < 4; k++) {
                total += m1.val[i][k] * m2.val[k][j];
            }

            result.val[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(const Matrix4 &m, Vec4 v) {
    double values[4];
    double total;

    for (int i = 0; i < 4; i++) {
        total = 0;
        for (int j = 0; j < 4; j++) {
            total += m.val[i][j] * v.getElementAt(j);
        }
        values[i] = total;
    }

    return {values[0], values[1], values[2], values[3], v.colorId};
}