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
#define PI 3.14159265

int culling(int i, Camera* cam, Triangle t, vector< vector<Vec3> > cameraview){
    
}

Matrix4 calculateViewportMatrix(Camera *camera) {
    Matrix4 res = getIdentityMatrix();

    res.val[0][0] = camera->verRes / 2.0;
    res.val[0][3] = ((camera->verRes - 1) / 2.0);
    res.val[1][1] = camera->horRes / 2.0;
    res.val[1][3] = ((camera->horRes - 1) / 2.0);
    res.val[2][2] = 0.5;
    res.val[2][3] = 0.5;

    res.val[3][3] = 0;

    return res;
}

bool LiangBarskyAlgorithm(Vec3 v0, Vec3 v1, Camera *cam, Vec3 &v0new, Vec3 &v1new) {
    double xmin = cam->left, xmax = cam->right, ymin = cam->bottom, ymax = cam->top, zmin = cam->near, zmax = cam->far;
    double te = 0.0, tl = 1.0;
    double dx = v1.x - v0.x, dy = v1.y - v0.y, dz = v1.z - v0.z;
    bool visible = false;
    v0new.x = v0.x;
    v0new.y = v0.y;
    v0new.z = v0.z;

    v1new.x = v1.x;
    v1new.y = v1.y;
    v1new.z = v1.z;

    if (isVisible(dx, xmin - v0.x, te, tl))
        if (isVisible(-dx, v0.x - xmax, te, tl))
            if (isVisible(dy, ymin - v0.y, te, tl))
                if (isVisible(-dy, v0.y - ymax, te, tl))
                    if (isVisible(dz, zmin - v0.z, te, tl))
                        if (isVisible(-dz, v0.z - zmax, te, tl)) {
                            visible = true;
                            if (tl < 1) {
                                v1new.x = v0.x + dx * tl;
                                v1new.y = v0.y + dy * tl;
                                v1new.z = v0.z + dz * tl;
                            }
                            if (te > 0) {
                                v0new.x = v0.x + dx * te;
                                v0new.y = v0.y + dy * te;
                                v0new.z = v0.z + dz * te;
                            }
                        }
    return visible;
}

bool isVisible(double den, double num, double &te, double &tl) {
    double t;
    if (den > 0) {
        t = num / den;
        if (t > tl) return false;
        if (t > te) te = t;
    } else if (den < 0) {
        t = num / den;
        if (t < te) return false;
        if (t < tl) tl = t;
    } else if (num > 0) return false; // parallel line
    return true;
}

Vec3 ComputeTriangleNormal(Vec3 v1, Vec3 v2, Vec3 v3) {
    Vec3 u = subtractVec3(v2, v1);
    Vec3 v = subtractVec3(v3, v1);

    return crossProductVec3(u, v);
}

bool IsTriangleVisible(Vec3 v1, Vec3 v2, Vec3 v3, Vec3 viewerPosition) {
    Vec3 normal = ComputeTriangleNormal(v1, v2, v3);
    Vec3 toViewer = subtractVec3(viewerPosition, v1);
// Check if a triangle is facing the viewer

    return dotProductVec3(normal, toViewer) > 0;
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

    //DONE: What if right-left = 0 -> Of course the plane would be a line if it happens!
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
    Vec3 u(r->ux, r->uy, r->uz, 0);
    Matrix4 fin;

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
    Vec3 w = crossProductVec3(u, v);
    double M[4][4] = {{u.x, u.y, u.z, 0},
                      {v.x, v.y, v.z, 0},
                      {w.x, w.y, w.z, 0},
                      {0,   0,   0,   1}};

    double inverse_M[4][4] = {{u.x, v.x, w.x, 0},
                              {u.y, v.y, w.y, 0},
                              {u.z, v.z, w.z, 0},
                              {0,   0,   0,   1}};

    double theta = rotation->angle * (PI / 180.0);
    double R_x[4][4] = {{1, 0,          0,           0},
                        {0, cos(theta), -sin(theta), 0},
                        {0, sin(theta), cos(theta),  0},
                        {0, 0,          0,           1}};


    Matrix4 r_x_m;
    r_x_m = multiplyMatrixWithMatrix(Matrix4(R_x), Matrix4(M));
    fin = multiplyMatrixWithMatrix(Matrix4(inverse_M), r_x_m);

    return fin;
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
Vec3 crossProductVec3(Vec3 a, Vec3 b) {
    Vec3 result;

    result.x = a.y * b.z - b.y * a.z;
    result.y = b.x * a.z - a.x * b.z;
    result.z = a.x * b.y - b.x * a.y;

    return result;
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v) {
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
Vec3 inverseVec3(Vec3 v) {
    Vec3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;

    return result;
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b) {
    Vec3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b) {
    Vec3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c) {
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
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2) {
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
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v) {
    double values[4];
    double total;

    for (int i = 0; i < 4; i++) {
        total = 0;
        for (int j = 0; j < 4; j++) {
            total += m.val[i][j] * v.getElementAt(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}