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


Matrix4 computeTranslationMatrix(Translation *t) {
    Matrix4 m = getIdentityMatrix();
    m.val[0][3] = t->tx;
    m.val[1][3] = t->ty;
    m.val[2][3] = t->tz;

    return m;
}
Matrix4 computeRotationMatrix(Rotation* rotation) {
        Rotation *r = rotation;
        double angle = r->angle;
        Vec3 u(r->ux, r->uy, r->uz, 0), v, w;
        Matrix4 fin;

        if (abs(u.x) <= abs(u.y) && abs(u.x) <= abs(u.z)) {
            v.x = 0;
            v.y = -u.z;
            v.z = u.y;
        }
        if (abs(u.y) <= abs(u.x) && abs(u.y) <= abs(u.z)) {
            v.x = -u.z;
            v.y = 0;
            v.z = u.x;
        }
        if (abs(u.x) <= abs(u.y) && abs(u.x) <= abs(u.z)) {
            v.x = -u.y;
            v.y = u.x;
            v.z = 0;
        }
        w = crossProductVec3(u, v);

        Matrix4 m = getIdentityMatrix(), mminus = getIdentityMatrix(), rx = getIdentityMatrix();
        m.val[0][0] = u.x;
        m.val[0][1] = u.y;
        m.val[0][2] = u.z;
        m.val[1][0] = v.x;
        m.val[1][1] = v.y;
        m.val[1][2] = v.z;
        m.val[2][0] = w.x;
        m.val[2][1] = w.y;
        m.val[2][2] = w.z;

        mminus.val[0][0] = u.x;
        mminus.val[0][1] = v.x;
        mminus.val[0][2] = w.x;
        mminus.val[0][0] = u.y;
        mminus.val[0][1] = v.x;
        mminus.val[0][2] = w.x;
        mminus.val[0][0] = u.z;
        mminus.val[0][1] = v.x;
        mminus.val[0][2] = w.x;

        rx.val[1][1] = cos(angle);
        rx.val[1][2] = -sin(angle);
        rx.val[2][1] = sin(angle);
        rx.val[2][2] = cos(angle);

        fin = multiplyMatrixWithMatrix(mminus, multiplyMatrixWithMatrix(rx, m));
        return fin;
}

Matrix4 computeScalingMatrix(Scaling *s) {
    Matrix4 m = getIdentityMatrix();
    m.val[0][3] = s->sx;
    m.val[1][3] = s->sy;
    m.val[2][3] = s->sz;

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