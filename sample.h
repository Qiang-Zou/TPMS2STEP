#ifndef SAMPLE_H
#define SAMPLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <Eigen/Core>

using namespace std;

void sample(std::vector<float>& params);
void sample2(std::vector<float>& params);
void sampleDiamond(std::vector<float>& params);
void sampleSchwarzP(std::vector<float>& params);
int findSpan2(int degree, int n, float t, float* knot);


#endif