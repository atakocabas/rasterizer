#ifndef __HELPERS_H__
#define __HELPERS_H__

#define ABS(a) ((a) > 0 ? (a) : -1 * (a))
#define EPSILON 0.000000001

#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Translation.h"
#include "Scaling.h"
#include "Rotation.h"
#include "Triangle.h"
#include "Camera.h"
#include "Color.h"
#include "Scene.h"
enum ORDER {
    BACK,
    FRONT
};

ORDER backfaceCulling(int i, int j, Camera *cam, vector< vector<Vec3> > cameraview);

bool LiangBarskyAlgorithm(const Vec3& v0, const Vec3& v1, Camera* cam, Vec3 &v0new, Vec3 &v1new, Color &color_a, Color &color_b, const Scene& scene);

bool isVisible(double den, double num, double &te, double &tl);

Matrix4 calculateViewportMatrix(Camera*);

Matrix4 calculatePerspectiveProjection(Camera* camera);

Matrix4 calculateOrthographicProjection(Camera *camera);
/*
WCS to CCS matrix.
*/
Matrix4 computeViewingTransformationMatrix(Camera*);
/*
* Compute all translation matrices.
*/
Matrix4 computeTranslationMatrix(Translation*);
/*
* Compute all rotation matrices.
*/
Matrix4 computeRotationMatrix(Rotation*);
/*
* Compute all scaling matrices.
*/
Matrix4 computeScalingMatrix(Scaling*);


/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(const Vec3& a, const Vec3& b);

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(const Vec3& a, const Vec3& b);

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(const Vec3& v);

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(const Vec3& v);

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(const Vec3& v);

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(const Vec3& a, const Vec3& b);

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(const Vec3& a, const Vec3& b);

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(const Vec3& v, double c);

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v);

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b);

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
*/
Matrix4 getIdentityMatrix();

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(const Matrix4& m1, const Matrix4& m2);

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(const Matrix4& m, Vec4 v);

#endif