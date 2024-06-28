#ifndef PIA_CUDA_H
#define PIA_CUDA_H

// #include <cuda.h>
#include "cuda_runtime.h"

#include <iostream>
#include <vector>
#include <ctime>
#include <chrono>
#include <cstdlib>
#include <cmath>

#include "Eigen/Core"

using namespace std;
using namespace std::chrono;
#define eps 1e-7

extern "C" void constrainedPIAGPU(float* knot, int degree, int numberOfPointCloud, int numberControlPointUDirection,
int numberControlPointVDirection, Eigen::MatrixXf& surfacePlusMatrix, Eigen::MatrixXf& surfaceMinusMatrix,
Eigen::VectorXf& U, Eigen::VectorXf& V, Eigen::MatrixXf* controlPointsPlusPIA, Eigen::MatrixXf* controlPointsMinusPIA,
float offsetValue, string modelType);

extern "C" void myTranslation(float* ctrlpts, int numberOfPoints, float tx, float ty, float tz);
extern "C" void clearMatrix(float* mat, int numberOfPoints);
extern "C" void Multiply(float* m, float* t, int pointNumber);
extern "C" void Multiply2(float* dest, float* src, float* t, int pointNumber);

extern "C" void memoryAllocation1(float* a, float* b, float* c, float* d, int number);
extern "C" void memoryAllocation2(float* a, float* b, float* c, float* d, float* e, int number);
extern "C" void myFree1(float* a, float* b, float* c, float* d);
extern "C" void myFree2(float* a, float* b, float* c, float* d, float* e);

extern "C" void initialzeCtrlPts(int number, float* ctrlpts, Eigen::MatrixXf& surface);

extern "C" void gyroidProcess();

#endif