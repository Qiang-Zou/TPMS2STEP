#include "piaCuda.h"

extern "C"
__device__ int findSpan(int degree, int n, float t, float* knot) {
    if (t == knot[n+1]) {
        return n;
    }
    int low = degree;
    int high = n + 1;
    int mid;

    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (t < knot[mid]-eps)
            high = mid;
        else
            low = mid;
    }
    return low;
}

extern "C"
__device__ void computeBasis(int degree, int span, float param, float* knot, float* basis) {
    float left[4] = {0};
    float right[4] = {0};
    basis[0] = 1.0;
    for (int j = 1; j <= degree; ++j) {
        left[j] = param - knot[span + 1 - j];
        right[j] = knot[span + j] - param;
        float saved = 0.0;
        for (int r = 0; r < j; ++r) {
            float basisMultiplier = basis[r] / (right[r + 1] + left[j - r]);
            basis[r] = saved + right[r + 1] * basisMultiplier;
            saved = left[j - r] * basisMultiplier;
        }
        basis[j] = saved;
    }
}

extern "C" 
void myTranslation(float* ctrlpts, int numberOfPoints, float tx, float ty, float tz) {
    for (int i=0;i<numberOfPoints;i++) {
        ctrlpts[i*3] += tx;
        ctrlpts[i*3+1] += ty;
        ctrlpts[i*3+2] += tz;
    }
    return;
}

extern "C"
void clearMatrix(float* mat, int numberOfPoints) {
    for (int i=0;i<numberOfPoints;i++) {
        mat[i*3] = 0;
        mat[i*3+1] = 0;
        mat[i*3+2] = 0;
    }
    return;
}

// NURBS surface evaluation
extern "C"
__global__ void evaluationGPU(float* knot, float* u, float* v, float* controlPoints, float* result,
int numberOfPointCloud, int numberControlPointUDirection, int numberControlPointVDirection) {
    int id = blockDim.x*blockIdx.x+threadIdx.x;
    if (id>=numberOfPointCloud) {
        return;
    }
    float r1 = 0, r2 = 0, r3 = 0;

    int spanU = findSpan(3, numberControlPointUDirection-1, u[id], knot);
    int spanV = findSpan(3, numberControlPointVDirection-1, v[id], knot);
    
    float basisU[4] = {0};
    float basisV[4] = {0};
    computeBasis(3, spanU, u[id], knot, basisU);
    computeBasis(3, spanV, v[id], knot, basisV);
    int uind = spanU - 3;
    int vind = spanV - 3;
    
    float temp[4][3];
    for (int l=0;l<=3;l++) {
        temp[l][0] = 0.0;
        temp[l][1] = 0.0;
        temp[l][2] = 0.0;
        vind = spanV-3+l;
        for (int k=0;k<=3;k++) {
            temp[l][0] = temp[l][0] + basisU[k]*controlPoints[(vind*numberControlPointUDirection+uind+k)*3];
            temp[l][1] = temp[l][1] + basisU[k]*controlPoints[(vind*numberControlPointUDirection+uind+k)*3+1];
            temp[l][2] = temp[l][2] + basisU[k]*controlPoints[(vind*numberControlPointUDirection+uind+k)*3+2];
        }
    }
    for (int q=0;q<=3;q++) {
        r1 += basisV[q]*temp[q][0];
        r2 += basisV[q]*temp[q][1];
        r3 += basisV[q]*temp[q][2];
    }
    result[id*3] = r1;
    result[id*3+1] = r2;
    result[id*3+2] = r3;
}

extern "C"
void Multiply(float* m, float* t, int pointNumber) {
    for (int i=0;i<pointNumber;i++) {
        float x = t[0]*m[i*3]+t[1]*m[i*3+1]+t[2]*m[i*3+2]+t[3];
        float y = t[4]*m[i*3]+t[5]*m[i*3+1]+t[6]*m[i*3+2]+t[7];
        float z = t[8]*m[i*3]+t[9]*m[i*3+1]+t[10]*m[i*3+2]+t[11];
        m[i*3] = x;
        m[i*3+1] = y;
        m[i*3+2] = z;
    }
}

extern "C"
void Multiply2(float* dest, float* src, float* t, int pointNumber) {
    for (int i=0;i<pointNumber;i++) {
        float x = t[0]*src[i*3]+t[1]*src[i*3+1]+t[2]*src[i*3+2]+t[3];
        float y = t[4]*src[i*3]+t[5]*src[i*3+1]+t[6]*src[i*3+2]+t[7];
        float z = t[8]*src[i*3]+t[9]*src[i*3+1]+t[10]*src[i*3+2]+t[11];
        dest[i*3] = x;
        dest[i*3+1] = y;
        dest[i*3+2] = z;
    }
}

extern "C"
void constrainedPIAGPU(float* knot, int degree, int numberOfPointCloud, int numberControlPointUDirection,
int numberControlPointVDirection, Eigen::MatrixXf& surfacePlusMatrix, Eigen::MatrixXf& surfaceMinusMatrix,
Eigen::VectorXf& U, Eigen::VectorXf& V, Eigen::MatrixXf* controlPointsPlusPIA, Eigen::MatrixXf* controlPointsMinusPIA,
float offsetValue, string modelType) {
    Eigen::MatrixXf evaluationSurfacePlusPIA(numberOfPointCloud, 3);
    Eigen::MatrixXf evaluationSurfaceMinusPIA(numberOfPointCloud, 3);

    // define variables for NURBS surface evaluation
    float* result = (float*)malloc(numberOfPointCloud*3*sizeof(float));
    float* uVector = (float*)malloc(numberOfPointCloud*sizeof(float));
    float* vVector = (float*)malloc(numberOfPointCloud*sizeof(float));

    for (int i=0;i<numberOfPointCloud;i++) {
        uVector[i] = U(i);
        vVector[i] = V(i);
    }

    int iterate = 20;
    int constraintIterationNumber = 6;

    // define the rigid transformation matrices
    float T4g[16], T5g[16], T6g[16], T7g[16];
    float T1d[16], T2d[16], T3d[16], T4d[16];
    float T1p[16], T2p[16], T3p[16], T4p[16];

    T4g[0] = 0;T4g[1] = 0;T4g[2] = 1;T4g[3] = 0.5;
    T4g[4] = 0;T4g[5] = -1;T4g[6] = 0;T4g[7] = 1.5;
    T4g[8] = 1;T4g[9] = 0;T4g[10] = 0;T4g[11] = -0.5;
    T4g[12] = 0;T4g[13] = 0;T4g[14] = 0;T4g[15] = 1;

    T5g[0] = -1;T5g[1] = 0;T5g[2] = 0;T5g[3] = 1.5;
    T5g[4] = 0;T5g[5] = 0;T5g[6] = -1;T5g[7] = 0.5;
    T5g[8] = 0;T5g[9] = -1;T5g[10] = 0;T5g[11] = 0.5;
    T5g[12] = 0;T5g[13] = 0;T5g[14] = 0;T5g[15] = 1;

    T6g[0] = 0;T6g[1] = 0;T6g[2] = -1;T6g[3] = 1.5;
    T6g[4] = 0;T6g[5] = -1;T6g[6] = 0;T6g[7] = 0.5;
    T6g[8] = -1;T6g[9] = 0;T6g[10] = 0;T6g[11] = 1.5;
    T6g[12] = 0;T6g[13] = 0;T6g[14] = 0;T6g[15] = 1;

    T7g[0] = -1;T7g[1] = 0;T7g[2] = 0;T7g[3] = 2.5;
    T7g[4] = 0;T7g[5] = 0;T7g[6] = 1;T7g[7] = 0.5;
    T7g[8] = 0;T7g[9] = 1;T7g[10] = 0;T7g[11] = -0.5;
    T7g[12] = 0;T7g[13] = 0;T7g[14] = 0;T7g[15] = 1;

    T1d[0] = 0;T1d[1] = -1;T1d[2] = 0;T1d[3] = 0;
    T1d[4] = 1;T1d[5] = 0;T1d[6] = 0;T1d[7] = 0;
    T1d[8] = 0;T1d[9] = 0;T1d[10] = -1;T1d[11] = 0;
    T1d[12] = 0;T1d[13] = 0;T1d[14] = 0;T1d[15] = 1;
    
    T2d[0] = -0.5;T2d[1] = -0.5;T2d[2] = sqrt(2)/2.0;T2d[3] = sqrt(2)/2.0;
    T2d[4] = -0.5;T2d[5] = -0.5;T2d[6] = -sqrt(2)/2.0;T2d[7] = sqrt(2)/2.0;
    T2d[8] = sqrt(2)/2.0;T2d[9] = -sqrt(2)/2.0;T2d[10] = 0;T2d[11] = 0;
    T2d[12] = 0;T2d[13] = 0;T2d[14] = 0;T2d[15] = 1;

    T3d[0] = -0.5;T3d[1] = 0.5;T3d[2] = sqrt(2)/2.0;T3d[3] = sqrt(2)/2.0;
    T3d[4] = 0.5;T3d[5] = -0.5;T3d[6] = sqrt(2)/2.0;T3d[7] = -sqrt(2)/2.0;
    T3d[8] = sqrt(2)/2.0;T3d[9] = sqrt(2)/2.0;T3d[10] = 0;T3d[11] = 0;
    T3d[12] = 0;T3d[13] = 0;T3d[14] = 0;T3d[15] = 1;

    T4d[0] = 0;T4d[1] = 1;T4d[2] = 0;T4d[3] = 0;
    T4d[4] = -1;T4d[5] = 0;T4d[6] = 0;T4d[7] = 0;
    T4d[8] = 0;T4d[9] = 0;T4d[10] = -1;T4d[11] = 0;
    T4d[12] = 0;T4d[13] = 0;T4d[14] = 0;T4d[15] = 1;

    T1p[0] = 0.5;T1p[1] = -0.5;T1p[2] = -sqrt(2)/2.0;T1p[3] = -sqrt(2)/4.0;
    T1p[4] = -0.5;T1p[5] = 0.5;T1p[6] = -sqrt(2)/2.0;T1p[7] = -sqrt(2)/4.0;
    T1p[8] = -sqrt(2)/2.0;T1p[9] = -sqrt(2)/2.0;T1p[10] = 0;T1p[11] = -0.5;
    T1p[12] = 0;T1p[13] = 0;T1p[14] = 0;T1p[15] = 1;

    T2p[0] = 0.5;T2p[1] = 0.5;T2p[2] = -sqrt(2)/2.0;T2p[3] = sqrt(2)/4.0;
    T2p[4] = 0.5;T2p[5] = 0.5;T2p[6] = sqrt(2)/2.0;T2p[7] = -sqrt(2)/4.0;
    T2p[8] = -sqrt(2)/2.0;T2p[9] = sqrt(2)/2.0;T2p[10] = 0;T2p[11] = 0.5;
    T2p[12] = 0;T2p[13] = 0;T2p[14] = 0;T2p[15] = 1;

    T3p[0] = 0;T3p[1] = -1;T3p[2] = 0;T3p[3] = 0;
    T3p[4] = -1;T3p[5] = 0;T3p[6] = 0;T3p[7] = 0;
    T3p[8] = 0;T3p[9] = 0;T3p[10] = 1;T3p[11] = 0;
    T3p[12] = 0;T3p[13] = 0;T3p[14] = 0;T3p[15] = 1;

    T4p[0] = 0;T4p[1] = 1;T4p[2] = 0;T4p[3] = 0;
    T4p[4] = 1;T4p[5] = 0;T4p[6] = 0;T4p[7] = 0;
    T4p[8] = 0;T4p[9] = 0;T4p[10] = 1;T4p[11] = 0;
    T4p[12] = 0;T4p[13] = 0;T4p[14] = 0;T4p[15] = 1;
    
    float* firstOrderEdge1 = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge1_origin = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge1_second = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge1_second_origin = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge1_minus = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge1_origin_minus = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge1_second_minus = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge1_second_origin_minus = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge2 = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge2_origin = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge2_second = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge2_second_origin = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge2_minus = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge2_origin_minus = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge2_second_minus = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge2_second_origin_minus = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge3 = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge3_origin = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge3_second = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge3_second_origin = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge3_minus = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge3_origin_minus = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge3_second_minus = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge3_second_origin_minus = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge4 = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge4_origin = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge4_second = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge4_second_origin = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge4_minus = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge4_origin_minus = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* firstOrderEdge4_second_minus = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
    float* firstOrderEdge4_second_origin_minus = (float*)malloc((numberControlPointVDirection-4)*3*sizeof(float));
    float* secondOrderEdge1 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_2 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_3 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_origin = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_second_origin = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_2_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_3_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_origin_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge1_second_origin_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_2 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_3 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_origin = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_second_origin = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_2_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_3_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_origin_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge2_second_origin_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_2 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_3 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_origin = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_second_origin = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_2_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_3_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_origin_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge3_second_origin_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_2 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_3 = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_origin = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_second_origin = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_2_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_3_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_origin_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
    float* secondOrderEdge4_second_origin_minus = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));

    // control points for two primitive surfaces with opposite offset
    float* controlPointsPlusPIAGPU = (float*)malloc(numberControlPointUDirection*numberControlPointVDirection*3*sizeof(float));
    float* controlPointsMinusPIAGPU = (float*)malloc(numberControlPointUDirection*numberControlPointVDirection*3*sizeof(float));
    
    // constraint
    float* constrainMatrix = (float*)malloc(numberControlPointUDirection*numberControlPointVDirection*3*sizeof(float));
    
    // initialize -- compute P0
    initialzeCtrlPts(numberControlPointUDirection*numberControlPointVDirection, controlPointsPlusPIAGPU, surfacePlusMatrix);
    initialzeCtrlPts(numberControlPointUDirection*numberControlPointVDirection, controlPointsMinusPIAGPU, surfaceMinusMatrix);

    dim3 threadsPerBlock(32);
    dim3 blocksPerGrid(std::ceil(numberOfPointCloud*1.0/threadsPerBlock.x));

    float* d_knot;
    float* d_u;
    float* d_v;
    float* d_controlPoints;
    float* d_result;

    cudaError_t cudastatus;
    cudastatus = cudaMalloc((void**)&d_knot, (numberControlPointUDirection+degree+1)*sizeof(float));
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > knot allocate error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    cudastatus = cudaMalloc((void**)&d_u, numberOfPointCloud*sizeof(float));
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > u allocate error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    cudastatus = cudaMalloc((void**)&d_v, numberOfPointCloud*sizeof(float));
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > v allocate error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    cudastatus = cudaMalloc((void**)&d_controlPoints, numberControlPointUDirection*numberControlPointVDirection*3*sizeof(float));
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > control points allocate error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    cudastatus = cudaMalloc((void**)&d_result, numberOfPointCloud*3*sizeof(float));
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > result allocate error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    cudastatus = cudaMemcpy(d_knot, knot, (numberControlPointUDirection+degree+1)*sizeof(float), cudaMemcpyHostToDevice);
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > knot transfer H2D error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    cudastatus = cudaMemcpy(d_u, uVector, numberOfPointCloud*sizeof(float), cudaMemcpyHostToDevice);
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > u transfer H2D error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    cudastatus = cudaMemcpy(d_v, vVector, numberOfPointCloud*sizeof(float), cudaMemcpyHostToDevice);
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > v transfer H2D error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    cudastatus = cudaMemcpy(d_controlPoints, controlPointsPlusPIAGPU, numberControlPointUDirection*numberControlPointVDirection*3*sizeof(float), cudaMemcpyHostToDevice);
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > control points 1 transfer H2D error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    // compute C
    evaluationGPU<<<blocksPerGrid, threadsPerBlock>>>(d_knot, d_u, d_v, d_controlPoints, d_result, numberOfPointCloud, numberControlPointUDirection, numberControlPointVDirection);

    cudastatus = cudaMemcpy(result, d_result, numberOfPointCloud*3*sizeof(float), cudaMemcpyDeviceToHost);
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > result 1 transfer D2H error: " << cudaGetErrorString(cudastatus) << std::endl;
    }
    
    for (int i=0;i<numberOfPointCloud;i++) {
        evaluationSurfacePlusPIA(i, 0) = result[i*3];
        evaluationSurfacePlusPIA(i, 1) = result[i*3+1];
        evaluationSurfacePlusPIA(i, 2) = result[i*3+2];
    }

    cudastatus = cudaMemcpy(d_controlPoints, controlPointsMinusPIAGPU, numberControlPointUDirection*numberControlPointVDirection*3*sizeof(float), cudaMemcpyHostToDevice);
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > control points 2 transfer H2D error: " << cudaGetErrorString(cudastatus) << std::endl;
    }

    evaluationGPU<<<blocksPerGrid, threadsPerBlock>>>(d_knot, d_u, d_v, d_controlPoints, d_result, numberOfPointCloud, numberControlPointUDirection, numberControlPointVDirection);
    
    cudastatus = cudaMemcpy(result, d_result, numberOfPointCloud*3*sizeof(float), cudaMemcpyDeviceToHost);
    if (cudaSuccess != cudastatus) {
        std::cout << "TPMS2STEP > result 2 transfer D2H error: " << cudaGetErrorString(cudastatus) << std::endl;
    }

    for (int i=0;i<numberOfPointCloud;i++) {
        evaluationSurfaceMinusPIA(i, 0) = result[i*3];
        evaluationSurfaceMinusPIA(i, 1) = result[i*3+1];
        evaluationSurfaceMinusPIA(i, 2) = result[i*3+2];
    }

    int count = 0;
    while (count < iterate) {
        count++;
        // compute P k+1
        clearMatrix(constrainMatrix, numberControlPointUDirection*numberControlPointVDirection);
        if (modelType == "Gyroid") {
            myTranslation(controlPointsPlusPIAGPU, numberControlPointUDirection*numberControlPointVDirection, 1, 0.5, 0);
            for (int i=0;i<numberControlPointUDirection;i++) {
                for (int j=0;j<numberControlPointVDirection;j++) {
                    if (j == numberControlPointVDirection-1 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        firstOrderEdge1[(i-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge1[(i-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge1[(i-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge1_second[(i-2)*3] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3];
                        firstOrderEdge1_second[(i-2)*3+1] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge1_second[(i-2)*3+2] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge1_origin[(i-2)*3] = firstOrderEdge1[(i-2)*3]-firstOrderEdge1_second[(i-2)*3];
                        firstOrderEdge1_origin[(i-2)*3+1] = firstOrderEdge1[(i-2)*3+1]-firstOrderEdge1_second[(i-2)*3+1];
                        firstOrderEdge1_origin[(i-2)*3+2] = firstOrderEdge1[(i-2)*3+2]-firstOrderEdge1_second[(i-2)*3+2];
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            secondOrderEdge1[(i-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge1[(i-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge1[(i-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge1_2[(i-3)*3] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3];
                            secondOrderEdge1_2[(i-3)*3+1] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge1_2[(i-3)*3+2] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge1_3[(i-3)*3] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3];
                            secondOrderEdge1_3[(i-3)*3+1] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge1_3[(i-3)*3+2] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge1_origin[(i-3)*3] = 2*secondOrderEdge1[(i-3)*3]-3*secondOrderEdge1_2[(i-3)*3]+secondOrderEdge1_3[(i-3)*3];
                            secondOrderEdge1_origin[(i-3)*3+1] = 2*secondOrderEdge1[(i-3)*3+1]-3*secondOrderEdge1_2[(i-3)*3+1]+secondOrderEdge1_3[(i-3)*3+1];
                            secondOrderEdge1_origin[(i-3)*3+2] = 2*secondOrderEdge1[(i-3)*3+2]-3*secondOrderEdge1_2[(i-3)*3+2]+secondOrderEdge1_3[(i-3)*3+2];
                        }
                    }
                    if (i == numberControlPointUDirection-1 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        firstOrderEdge2[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge2[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge2[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge2_second[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3];
                        firstOrderEdge2_second[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+1];
                        firstOrderEdge2_second[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+2];
                        firstOrderEdge2_origin[(j-2)*3] = firstOrderEdge2[(j-2)*3]-firstOrderEdge2_second[(j-2)*3];
                        firstOrderEdge2_origin[(j-2)*3+1] = firstOrderEdge2[(j-2)*3+1]-firstOrderEdge2_second[(j-2)*3+1];
                        firstOrderEdge2_origin[(j-2)*3+2] = firstOrderEdge2[(j-2)*3+2]-firstOrderEdge2_second[(j-2)*3+2];
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            secondOrderEdge2[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge2[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_2[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3];
                            secondOrderEdge2_2[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+1];
                            secondOrderEdge2_2[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+2];
                            secondOrderEdge2_3[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3];
                            secondOrderEdge2_3[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3+1];
                            secondOrderEdge2_3[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3+2];
                            secondOrderEdge2_origin[(j-3)*3] = 2*secondOrderEdge2[(j-3)*3]-3*secondOrderEdge2_2[(j-3)*3]+secondOrderEdge2_3[(j-3)*3];
                            secondOrderEdge2_origin[(j-3)*3+1] = 2*secondOrderEdge2[(j-3)*3+1]-3*secondOrderEdge2_2[(j-3)*3+1]+secondOrderEdge2_3[(j-3)*3+1];
                            secondOrderEdge2_origin[(j-3)*3+2] = 2*secondOrderEdge2[(j-3)*3+2]-3*secondOrderEdge2_2[(j-3)*3+2]+secondOrderEdge2_3[(j-3)*3+2];
                        }
                    }
                    if (j == 0 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        firstOrderEdge3[(i-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge3[(i-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge3[(i-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge3_second[(i-2)*3] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3];
                        firstOrderEdge3_second[(i-2)*3+1] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge3_second[(i-2)*3+2] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge3_origin[(i-2)*3] = firstOrderEdge3[(i-2)*3]-firstOrderEdge3_second[(i-2)*3];
                        firstOrderEdge3_origin[(i-2)*3+1] = firstOrderEdge3[(i-2)*3+1]-firstOrderEdge3_second[(i-2)*3+1];
                        firstOrderEdge3_origin[(i-2)*3+2] = firstOrderEdge3[(i-2)*3+2]-firstOrderEdge3_second[(i-2)*3+2];
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            secondOrderEdge3[(i-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge3[(i-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge3[(i-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge3_2[(i-3)*3] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3];
                            secondOrderEdge3_2[(i-3)*3+1] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge3_2[(i-3)*3+2] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge3_3[(i-3)*3] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3];
                            secondOrderEdge3_3[(i-3)*3+1] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge3_3[(i-3)*3+2] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge3_origin[(i-3)*3] = 2*secondOrderEdge3[(i-3)*3]-3*secondOrderEdge3_2[(i-3)*3]+secondOrderEdge3_3[(i-3)*3];
                            secondOrderEdge3_origin[(i-3)*3+1] = 2*secondOrderEdge3[(i-3)*3+1]-3*secondOrderEdge3_2[(i-3)*3+1]+secondOrderEdge3_3[(i-3)*3+1];
                            secondOrderEdge3_origin[(i-3)*3+2] = 2*secondOrderEdge3[(i-3)*3+2]-3*secondOrderEdge3_2[(i-3)*3+2]+secondOrderEdge3_3[(i-3)*3+2];
                        }
                    }
                    if (i == 0 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        firstOrderEdge4[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge4[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge4[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];

                        firstOrderEdge4_second[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3];
                        firstOrderEdge4_second[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+1];
                        firstOrderEdge4_second[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+2];

                        firstOrderEdge4_origin[(j-2)*3] = firstOrderEdge4[(j-2)*3]-firstOrderEdge4_second[(j-2)*3];
                        firstOrderEdge4_origin[(j-2)*3+1] = firstOrderEdge4[(j-2)*3+1]-firstOrderEdge4_second[(j-2)*3+1];
                        firstOrderEdge4_origin[(j-2)*3+2] = firstOrderEdge4[(j-2)*3+2]-firstOrderEdge4_second[(j-2)*3+2];

                        if (j != 2 && j != numberControlPointVDirection-3) {
                            secondOrderEdge4[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge4[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_2[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3];
                            secondOrderEdge4_2[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+1];
                            secondOrderEdge4_2[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+2];
                            secondOrderEdge4_3[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3];
                            secondOrderEdge4_3[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3+1];
                            secondOrderEdge4_3[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3+2];
                            secondOrderEdge4_origin[(j-3)*3] = 2*secondOrderEdge4[(j-3)*3]-3*secondOrderEdge4_2[(j-3)*3]+secondOrderEdge4_3[(j-3)*3];
                            secondOrderEdge4_origin[(j-3)*3+1] = 2*secondOrderEdge4[(j-3)*3+1]-3*secondOrderEdge4_2[(j-3)*3+1]+secondOrderEdge4_3[(j-3)*3+1];
                            secondOrderEdge4_origin[(j-3)*3+2] = 2*secondOrderEdge4[(j-3)*3+2]-3*secondOrderEdge4_2[(j-3)*3+2]+secondOrderEdge4_3[(j-3)*3+2];
                        }
                    }
                }
            }

            // for (int i=0;i<numberControlPointUDirection-4;i++) {
            //     cout << "1: " << firstOrderEdge1[i*3] << " " << firstOrderEdge1[i*3+1] << " " << firstOrderEdge1[i*3+2] << endl;
            // }

            // rigid transformation
            Multiply(firstOrderEdge1, T4g, numberControlPointUDirection-4);
            Multiply(firstOrderEdge1_second, T4g, numberControlPointUDirection-4);
            Multiply(secondOrderEdge1, T4g, numberControlPointUDirection-6);
            Multiply(secondOrderEdge1_2, T4g, numberControlPointUDirection-6);
            Multiply(secondOrderEdge1_3, T4g, numberControlPointUDirection-6);
            Multiply(firstOrderEdge2, T5g, numberControlPointUDirection-4);
            Multiply(firstOrderEdge2_second, T5g, numberControlPointUDirection-4);
            Multiply(secondOrderEdge2, T5g, numberControlPointUDirection-6);
            Multiply(secondOrderEdge2_2, T5g, numberControlPointUDirection-6);
            Multiply(secondOrderEdge2_3, T5g, numberControlPointUDirection-6);
            Multiply(firstOrderEdge3, T6g, numberControlPointUDirection-4);
            Multiply(firstOrderEdge3_second, T6g, numberControlPointUDirection-4);
            Multiply(secondOrderEdge3, T6g, numberControlPointUDirection-6);
            Multiply(secondOrderEdge3_2, T6g, numberControlPointUDirection-6);
            Multiply(secondOrderEdge3_3, T6g, numberControlPointUDirection-6);
            Multiply(firstOrderEdge4, T7g, numberControlPointUDirection-4);
            Multiply(firstOrderEdge4_second, T7g, numberControlPointUDirection-4);
            Multiply(secondOrderEdge4, T7g, numberControlPointUDirection-6);
            Multiply(secondOrderEdge4_2, T7g, numberControlPointUDirection-6);
            Multiply(secondOrderEdge4_3, T7g, numberControlPointUDirection-6);

            // for (int i=0;i<numberControlPointUDirection-4;i++) {
            //     cout << "2: " << firstOrderEdge1[i*3] << " " << firstOrderEdge1[i*3+1] << " " << firstOrderEdge1[i*3+2] << endl;
            // }

            for (int i=0;i<numberControlPointUDirection-4;i++) {
                firstOrderEdge1_second_origin[i*3] = (-1)*(firstOrderEdge1[i*3]-firstOrderEdge1_second[i*3]);
                firstOrderEdge1_second_origin[i*3+1] = (-1)*(firstOrderEdge1[i*3+1]-firstOrderEdge1_second[i*3+1]);
                firstOrderEdge1_second_origin[i*3+2] = (-1)*(firstOrderEdge1[i*3+2]-firstOrderEdge1_second[i*3+2]);
                firstOrderEdge2_second_origin[i*3] = (-1)*(firstOrderEdge2[i*3]-firstOrderEdge2_second[i*3]);
                firstOrderEdge2_second_origin[i*3+1] = (-1)*(firstOrderEdge2[i*3+1]-firstOrderEdge2_second[i*3+1]);
                firstOrderEdge2_second_origin[i*3+2] = (-1)*(firstOrderEdge2[i*3+2]-firstOrderEdge2_second[i*3+2]);
                firstOrderEdge3_second_origin[i*3] = (-1)*(firstOrderEdge3[i*3]-firstOrderEdge3_second[i*3]);
                firstOrderEdge3_second_origin[i*3+1] = (-1)*(firstOrderEdge3[i*3+1]-firstOrderEdge3_second[i*3+1]);
                firstOrderEdge3_second_origin[i*3+2] = (-1)*(firstOrderEdge3[i*3+2]-firstOrderEdge3_second[i*3+2]);
                firstOrderEdge4_second_origin[i*3] = (-1)*(firstOrderEdge4[i*3]-firstOrderEdge4_second[i*3]);
                firstOrderEdge4_second_origin[i*3+1] = (-1)*(firstOrderEdge4[i*3+1]-firstOrderEdge4_second[i*3+1]);
                firstOrderEdge4_second_origin[i*3+2] = (-1)*(firstOrderEdge4[i*3+2]-firstOrderEdge4_second[i*3+2]);
            }
            for (int i=0;i<numberControlPointUDirection-6;i++) {
                secondOrderEdge1_second_origin[i*3] = (-1)*(2*secondOrderEdge1[i*3]-3*secondOrderEdge1_2[i*3]+secondOrderEdge1_3[i*3]);
                secondOrderEdge1_second_origin[i*3+1] = (-1)*(2*secondOrderEdge1[i*3+1]-3*secondOrderEdge1_2[i*3+1]+secondOrderEdge1_3[i*3+1]);
                secondOrderEdge1_second_origin[i*3+2] = (-1)*(2*secondOrderEdge1[i*3+2]-3*secondOrderEdge1_2[i*3+2]+secondOrderEdge1_3[i*3+2]);
                secondOrderEdge2_second_origin[i*3] = (-1)*(2*secondOrderEdge2[i*3]-3*secondOrderEdge2_2[i*3]+secondOrderEdge2_3[i*3]);
                secondOrderEdge2_second_origin[i*3+1] = (-1)*(2*secondOrderEdge2[i*3+1]-3*secondOrderEdge2_2[i*3+1]+secondOrderEdge2_3[i*3+1]);
                secondOrderEdge2_second_origin[i*3+2] = (-1)*(2*secondOrderEdge2[i*3+2]-3*secondOrderEdge2_2[i*3+2]+secondOrderEdge2_3[i*3+2]);
                secondOrderEdge3_second_origin[i*3] = (-1)*(2*secondOrderEdge3[i*3]-3*secondOrderEdge3_2[i*3]+secondOrderEdge3_3[i*3]);
                secondOrderEdge3_second_origin[i*3+1] = (-1)*(2*secondOrderEdge3[i*3+1]-3*secondOrderEdge3_2[i*3+1]+secondOrderEdge3_3[i*3+1]);
                secondOrderEdge3_second_origin[i*3+2] = (-1)*(2*secondOrderEdge3[i*3+2]-3*secondOrderEdge3_2[i*3+2]+secondOrderEdge3_3[i*3+2]);
                secondOrderEdge4_second_origin[i*3] = (-1)*(2*secondOrderEdge4[i*3]-3*secondOrderEdge4_2[i*3]+secondOrderEdge4_3[i*3]);
                secondOrderEdge4_second_origin[i*3+1] = (-1)*(2*secondOrderEdge4[i*3+1]-3*secondOrderEdge4_2[i*3+1]+secondOrderEdge4_3[i*3+1]);
                secondOrderEdge4_second_origin[i*3+2] = (-1)*(2*secondOrderEdge4[i*3+2]-3*secondOrderEdge4_2[i*3+2]+secondOrderEdge4_3[i*3+2]);
            }

            for (int i=0;i<numberControlPointUDirection;i++) {
                for (int j=0;j<numberControlPointVDirection;j++) {
                    if (j == numberControlPointVDirection-2 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge1_second_origin[(numberControlPointUDirection-5-(i-2))*3]-firstOrderEdge1_origin[(i-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge1_second_origin[(numberControlPointUDirection-5-(i-2))*3+1]-firstOrderEdge1_origin[(i-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge1_second_origin[(numberControlPointUDirection-5-(i-2))*3+2]-firstOrderEdge1_origin[(i-2)*3+2])/2;
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3] = (secondOrderEdge1_second_origin[(numberControlPointUDirection-7-(i-3))*3]-secondOrderEdge1_origin[(i-3)*3])/12.0;
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3+1] = (secondOrderEdge1_second_origin[(numberControlPointUDirection-7-(i-3))*3+1]-secondOrderEdge1_origin[(i-3)*3+1])/12.0;
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3+2] = (secondOrderEdge1_second_origin[(numberControlPointUDirection-7-(i-3))*3+2]-secondOrderEdge1_origin[(i-3)*3+2])/12.0;
                        }
                    }
                    if (i == numberControlPointUDirection-2 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge2_second_origin[(numberControlPointUDirection-5-(j-2))*3]-firstOrderEdge2_origin[(j-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge2_second_origin[(numberControlPointUDirection-5-(j-2))*3+1]-firstOrderEdge2_origin[(j-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge2_second_origin[(numberControlPointUDirection-5-(j-2))*3+2]-firstOrderEdge2_origin[(j-2)*3+2])/2;
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3] = (secondOrderEdge2_second_origin[(numberControlPointUDirection-7-(j-3))*3]-secondOrderEdge2_origin[(j-3)*3])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3+1] = (secondOrderEdge2_second_origin[(numberControlPointUDirection-7-(j-3))*3+1]-secondOrderEdge2_origin[(j-3)*3+1])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3+2] = (secondOrderEdge2_second_origin[(numberControlPointUDirection-7-(j-3))*3+2]-secondOrderEdge2_origin[(j-3)*3+2])/12.0;
                        }
                    }
                    if (j == 1 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge3_second_origin[(numberControlPointUDirection-5-(i-2))*3]-firstOrderEdge3_origin[(i-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge3_second_origin[(numberControlPointUDirection-5-(i-2))*3+1]-firstOrderEdge3_origin[(i-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge3_second_origin[(numberControlPointUDirection-5-(i-2))*3+2]-firstOrderEdge3_origin[(i-2)*3+2])/2;
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3] = (secondOrderEdge3_second_origin[(numberControlPointUDirection-7-(i-3))*3]-secondOrderEdge3_origin[(i-3)*3])/12.0;
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3+1] = (secondOrderEdge3_second_origin[(numberControlPointUDirection-7-(i-3))*3+1]-secondOrderEdge3_origin[(i-3)*3+1])/12.0;
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3+2] = (secondOrderEdge3_second_origin[(numberControlPointUDirection-7-(i-3))*3+2]-secondOrderEdge3_origin[(i-3)*3+2])/12.0;
                        }
                    }
                    if (i == 1 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge4_second_origin[(numberControlPointUDirection-5-(j-2))*3]-firstOrderEdge4_origin[(j-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge4_second_origin[(numberControlPointUDirection-5-(j-2))*3+1]-firstOrderEdge4_origin[(j-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge4_second_origin[(numberControlPointUDirection-5-(j-2))*3+2]-firstOrderEdge4_origin[(j-2)*3+2])/2;
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3] = (secondOrderEdge4_second_origin[(numberControlPointUDirection-7-(j-3))*3]-secondOrderEdge4_origin[(j-3)*3])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3+1] = (secondOrderEdge4_second_origin[(numberControlPointUDirection-7-(j-3))*3+1]-secondOrderEdge4_origin[(j-3)*3+1])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3+2] = (secondOrderEdge4_second_origin[(numberControlPointUDirection-7-(j-3))*3+2]-secondOrderEdge4_origin[(j-3)*3+2])/12.0;
                        }
                    }
                }
            }
            myTranslation(controlPointsPlusPIAGPU, numberControlPointUDirection*numberControlPointVDirection, -1, -0.5, 0);
        }
        else if (modelType == "Diamond") {
            clearMatrix(constrainMatrix, numberControlPointUDirection*numberControlPointVDirection);
            for (int i=0;i<numberControlPointUDirection;i++) {
                for (int j=0;j<numberControlPointVDirection;j++) {
                    //edge 2
                    if (j == 0 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        firstOrderEdge2[(i-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge2[(i-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge2[(i-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge2_second[(i-2)*3] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3];
                        firstOrderEdge2_second[(i-2)*3+1] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge2_second[(i-2)*3+2] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge2_origin[(i-2)*3] = firstOrderEdge2[(i-2)*3]-firstOrderEdge2_second[(i-2)*3];
                        firstOrderEdge2_origin[(i-2)*3+1] = firstOrderEdge2[(i-2)*3+1]-firstOrderEdge2_second[(i-2)*3+1];
                        firstOrderEdge2_origin[(i-2)*3+2] = firstOrderEdge2[(i-2)*3+2]-firstOrderEdge2_second[(i-2)*3+2];
                        firstOrderEdge2_minus[(i-2)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge2_minus[(i-2)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge2_minus[(i-2)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge2_second_minus[(i-2)*3] = controlPointsMinusPIAGPU[((j+1)*numberControlPointUDirection+i)*3];
                        firstOrderEdge2_second_minus[(i-2)*3+1] = controlPointsMinusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge2_second_minus[(i-2)*3+2] = controlPointsMinusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge2_origin_minus[(i-2)*3] = firstOrderEdge2_minus[(i-2)*3]-firstOrderEdge2_second_minus[(i-2)*3];
                        firstOrderEdge2_origin_minus[(i-2)*3+1] = firstOrderEdge2_minus[(i-2)*3+1]-firstOrderEdge2_second_minus[(i-2)*3+1];
                        firstOrderEdge2_origin_minus[(i-2)*3+2] = firstOrderEdge2_minus[(i-2)*3+2]-firstOrderEdge2_second_minus[(i-2)*3+2];
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            secondOrderEdge2[(i-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge2[(i-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2[(i-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_2[(i-3)*3] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3];
                            secondOrderEdge2_2[(i-3)*3+1] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2_2[(i-3)*3+2] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_3[(i-3)*3] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3];
                            secondOrderEdge2_3[(i-3)*3+1] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2_3[(i-3)*3+2] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_origin[(i-3)*3] = 2*secondOrderEdge2[(i-3)*3]-3*secondOrderEdge2_2[(i-3)*3]+secondOrderEdge2_3[(i-3)*3];
                            secondOrderEdge2_origin[(i-3)*3+1] = 2*secondOrderEdge2[(i-3)*3+1]-3*secondOrderEdge2_2[(i-3)*3+1]+secondOrderEdge2_3[(i-3)*3+1];
                            secondOrderEdge2_origin[(i-3)*3+2] = 2*secondOrderEdge2[(i-3)*3+2]-3*secondOrderEdge2_2[(i-3)*3+2]+secondOrderEdge2_3[(i-3)*3+2];
                            secondOrderEdge2_minus[(i-3)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge2_minus[(i-3)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2_minus[(i-3)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_2_minus[(i-3)*3] = controlPointsMinusPIAGPU[((j+1)*numberControlPointUDirection+i)*3];
                            secondOrderEdge2_2_minus[(i-3)*3+1] = controlPointsMinusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2_2_minus[(i-3)*3+2] = controlPointsMinusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_3_minus[(i-3)*3] = controlPointsMinusPIAGPU[((j+2)*numberControlPointUDirection+i)*3];
                            secondOrderEdge2_3_minus[(i-3)*3+1] = controlPointsMinusPIAGPU[((j+2)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2_3_minus[(i-3)*3+2] = controlPointsMinusPIAGPU[((j+2)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_origin_minus[(i-3)*3] = 2*secondOrderEdge2_minus[(i-3)*3]-3*secondOrderEdge2_2_minus[(i-3)*3]+secondOrderEdge2_3_minus[(i-3)*3];
                            secondOrderEdge2_origin_minus[(i-3)*3+1] = 2*secondOrderEdge2_minus[(i-3)*3+1]-3*secondOrderEdge2_2_minus[(i-3)*3+1]+secondOrderEdge2_3_minus[(i-3)*3+1];
                            secondOrderEdge2_origin_minus[(i-3)*3+2] = 2*secondOrderEdge2_minus[(i-3)*3+2]-3*secondOrderEdge2_2_minus[(i-3)*3+2]+secondOrderEdge2_3_minus[(i-3)*3+2];
                        }
                    }
                    //edge 1
                    if (i == numberControlPointUDirection-1 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        firstOrderEdge1[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge1[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge1[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge1_second[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3];
                        firstOrderEdge1_second[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+1];
                        firstOrderEdge1_second[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+2];
                        firstOrderEdge1_origin[(j-2)*3] = firstOrderEdge1[(j-2)*3]-firstOrderEdge1_second[(j-2)*3];
                        firstOrderEdge1_origin[(j-2)*3+1] = firstOrderEdge1[(j-2)*3+1]-firstOrderEdge1_second[(j-2)*3+1];
                        firstOrderEdge1_origin[(j-2)*3+2] = firstOrderEdge1[(j-2)*3+2]-firstOrderEdge1_second[(j-2)*3+2];
                        firstOrderEdge1_minus[(j-2)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge1_minus[(j-2)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge1_minus[(j-2)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge1_second_minus[(j-2)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-1))*3];
                        firstOrderEdge1_second_minus[(j-2)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+1];
                        firstOrderEdge1_second_minus[(j-2)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+2];
                        firstOrderEdge1_origin_minus[(j-2)*3] = firstOrderEdge1_minus[(j-2)*3]-firstOrderEdge1_second_minus[(j-2)*3];
                        firstOrderEdge1_origin_minus[(j-2)*3+1] = firstOrderEdge1_minus[(j-2)*3+1]-firstOrderEdge1_second_minus[(j-2)*3+1];
                        firstOrderEdge1_origin_minus[(j-2)*3+2] = firstOrderEdge1_minus[(j-2)*3+2]-firstOrderEdge1_second_minus[(j-2)*3+2];
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            secondOrderEdge1[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge1[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge1[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge1_2[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3];
                            secondOrderEdge1_2[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+1];
                            secondOrderEdge1_2[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+2];
                            secondOrderEdge1_3[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3];
                            secondOrderEdge1_3[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3+1];
                            secondOrderEdge1_3[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3+2];
                            secondOrderEdge1_origin[(j-3)*3] = 2*secondOrderEdge1[(j-3)*3]-3*secondOrderEdge1_2[(j-3)*3]+secondOrderEdge1_3[(j-3)*3];
                            secondOrderEdge1_origin[(j-3)*3+1] = 2*secondOrderEdge1[(j-3)*3+1]-3*secondOrderEdge1_2[(j-3)*3+1]+secondOrderEdge1_3[(j-3)*3+1];
                            secondOrderEdge1_origin[(j-3)*3+2] = 2*secondOrderEdge1[(j-3)*3+2]-3*secondOrderEdge1_2[(j-3)*3+2]+secondOrderEdge1_3[(j-3)*3+2];
                            secondOrderEdge1_minus[(j-3)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge1_minus[(j-3)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge1_minus[(j-3)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge1_2_minus[(j-3)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-1))*3];
                            secondOrderEdge1_2_minus[(j-3)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+1];
                            secondOrderEdge1_2_minus[(j-3)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+2];
                            secondOrderEdge1_3_minus[(j-3)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-2))*3];
                            secondOrderEdge1_3_minus[(j-3)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-2))*3+1];
                            secondOrderEdge1_3_minus[(j-3)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i-2))*3+2];
                            secondOrderEdge1_origin_minus[(j-3)*3] = 2*secondOrderEdge1_minus[(j-3)*3]-3*secondOrderEdge1_2_minus[(j-3)*3]+secondOrderEdge1_3_minus[(j-3)*3];
                            secondOrderEdge1_origin_minus[(j-3)*3+1] = 2*secondOrderEdge1_minus[(j-3)*3+1]-3*secondOrderEdge1_2_minus[(j-3)*3+1]+secondOrderEdge1_3_minus[(j-3)*3+1];
                            secondOrderEdge1_origin_minus[(j-3)*3+2] = 2*secondOrderEdge1_minus[(j-3)*3+2]-3*secondOrderEdge1_2_minus[(j-3)*3+2]+secondOrderEdge1_3_minus[(j-3)*3+2];
                        }
                    }
                    // edge 4
                    if (j == numberControlPointVDirection-1 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        firstOrderEdge4[(i-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge4[(i-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge4[(i-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge4_second[(i-2)*3] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3];
                        firstOrderEdge4_second[(i-2)*3+1] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge4_second[(i-2)*3+2] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge4_origin[(i-2)*3] = firstOrderEdge4[(i-2)*3]-firstOrderEdge4_second[(i-2)*3];
                        firstOrderEdge4_origin[(i-2)*3+1] = firstOrderEdge4[(i-2)*3+1]-firstOrderEdge4_second[(i-2)*3+1];
                        firstOrderEdge4_origin[(i-2)*3+2] = firstOrderEdge4[(i-2)*3+2]-firstOrderEdge4_second[(i-2)*3+2];
                        firstOrderEdge4_minus[(i-2)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge4_minus[(i-2)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge4_minus[(i-2)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge4_second_minus[(i-2)*3] = controlPointsMinusPIAGPU[((j-1)*numberControlPointUDirection+i)*3];
                        firstOrderEdge4_second_minus[(i-2)*3+1] = controlPointsMinusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge4_second_minus[(i-2)*3+2] = controlPointsMinusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge4_origin_minus[(i-2)*3] = firstOrderEdge4_minus[(i-2)*3]-firstOrderEdge4_second_minus[(i-2)*3];
                        firstOrderEdge4_origin_minus[(i-2)*3+1] = firstOrderEdge4_minus[(i-2)*3+1]-firstOrderEdge4_second_minus[(i-2)*3+1];
                        firstOrderEdge4_origin_minus[(i-2)*3+2] = firstOrderEdge4_minus[(i-2)*3+2]-firstOrderEdge4_second_minus[(i-2)*3+2];
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            secondOrderEdge4[(i-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge4[(i-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4[(i-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_2[(i-3)*3] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3];
                            secondOrderEdge4_2[(i-3)*3+1] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4_2[(i-3)*3+2] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_3[(i-3)*3] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3];
                            secondOrderEdge4_3[(i-3)*3+1] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4_3[(i-3)*3+2] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_origin[(i-3)*3] = 2*secondOrderEdge4[(i-3)*3]-3*secondOrderEdge4_2[(i-3)*3]+secondOrderEdge4_3[(i-3)*3];
                            secondOrderEdge4_origin[(i-3)*3+1] = 2*secondOrderEdge4[(i-3)*3+1]-3*secondOrderEdge4_2[(i-3)*3+1]+secondOrderEdge4_3[(i-3)*3+1];
                            secondOrderEdge4_origin[(i-3)*3+2] = 2*secondOrderEdge4[(i-3)*3+2]-3*secondOrderEdge4_2[(i-3)*3+2]+secondOrderEdge4_3[(i-3)*3+2];
                            secondOrderEdge4_minus[(i-3)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge4_minus[(i-3)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4_minus[(i-3)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_2_minus[(i-3)*3] = controlPointsMinusPIAGPU[((j-1)*numberControlPointUDirection+i)*3];
                            secondOrderEdge4_2_minus[(i-3)*3+1] = controlPointsMinusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4_2_minus[(i-3)*3+2] = controlPointsMinusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_3_minus[(i-3)*3] = controlPointsMinusPIAGPU[((j-2)*numberControlPointUDirection+i)*3];
                            secondOrderEdge4_3_minus[(i-3)*3+1] = controlPointsMinusPIAGPU[((j-2)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4_3_minus[(i-3)*3+2] = controlPointsMinusPIAGPU[((j-2)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_origin_minus[(i-3)*3] = 2*secondOrderEdge4_minus[(i-3)*3]-3*secondOrderEdge4_2_minus[(i-3)*3]+secondOrderEdge4_3_minus[(i-3)*3];
                            secondOrderEdge4_origin_minus[(i-3)*3+1] = 2*secondOrderEdge4_minus[(i-3)*3+1]-3*secondOrderEdge4_2_minus[(i-3)*3+1]+secondOrderEdge4_3_minus[(i-3)*3+1];
                            secondOrderEdge4_origin_minus[(i-3)*3+2] = 2*secondOrderEdge4_minus[(i-3)*3+2]-3*secondOrderEdge4_2_minus[(i-3)*3+2]+secondOrderEdge4_3_minus[(i-3)*3+2];
                        }
                    }
                    // edge 3
                    if (i == 0 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        firstOrderEdge3[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge3[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge3[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge3_second[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3];
                        firstOrderEdge3_second[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+1];
                        firstOrderEdge3_second[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+2];
                        firstOrderEdge3_origin[(j-2)*3] = firstOrderEdge3[(j-2)*3]-firstOrderEdge3_second[(j-2)*3];
                        firstOrderEdge3_origin[(j-2)*3+1] = firstOrderEdge3[(j-2)*3+1]-firstOrderEdge3_second[(j-2)*3+1];
                        firstOrderEdge3_origin[(j-2)*3+2] = firstOrderEdge3[(j-2)*3+2]-firstOrderEdge3_second[(j-2)*3+2];
                        firstOrderEdge3_minus[(j-2)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge3_minus[(j-2)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge3_minus[(j-2)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge3_second_minus[(j-2)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+1))*3];
                        firstOrderEdge3_second_minus[(j-2)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+1];
                        firstOrderEdge3_second_minus[(j-2)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+2];
                        firstOrderEdge3_origin_minus[(j-2)*3] = firstOrderEdge3_minus[(j-2)*3]-firstOrderEdge3_second_minus[(j-2)*3];
                        firstOrderEdge3_origin_minus[(j-2)*3+1] = firstOrderEdge3_minus[(j-2)*3+1]-firstOrderEdge3_second_minus[(j-2)*3+1];
                        firstOrderEdge3_origin_minus[(j-2)*3+2] = firstOrderEdge3_minus[(j-2)*3+2]-firstOrderEdge3_second_minus[(j-2)*3+2];
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            secondOrderEdge3[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge3[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge3[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge3_2[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3];
                            secondOrderEdge3_2[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+1];
                            secondOrderEdge3_2[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+2];
                            secondOrderEdge3_3[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3];
                            secondOrderEdge3_3[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3+1];
                            secondOrderEdge3_3[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3+2];
                            secondOrderEdge3_origin[(j-3)*3] = 2*secondOrderEdge3[(j-3)*3]-3*secondOrderEdge3_2[(j-3)*3]+secondOrderEdge3_3[(j-3)*3];
                            secondOrderEdge3_origin[(j-3)*3+1] = 2*secondOrderEdge3[(j-3)*3+1]-3*secondOrderEdge3_2[(j-3)*3+1]+secondOrderEdge3_3[(j-3)*3+1];
                            secondOrderEdge3_origin[(j-3)*3+2] = 2*secondOrderEdge3[(j-3)*3+2]-3*secondOrderEdge3_2[(j-3)*3+2]+secondOrderEdge3_3[(j-3)*3+2];
                            secondOrderEdge3_minus[(j-3)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge3_minus[(j-3)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge3_minus[(j-3)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge3_2_minus[(j-3)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+1))*3];
                            secondOrderEdge3_2_minus[(j-3)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+1];
                            secondOrderEdge3_2_minus[(j-3)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+2];
                            secondOrderEdge3_3_minus[(j-3)*3] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+2))*3];
                            secondOrderEdge3_3_minus[(j-3)*3+1] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+2))*3+1];
                            secondOrderEdge3_3_minus[(j-3)*3+2] = controlPointsMinusPIAGPU[(j*numberControlPointUDirection+(i+2))*3+2];
                            secondOrderEdge3_origin_minus[(j-3)*3] = 2*secondOrderEdge3_minus[(j-3)*3]-3*secondOrderEdge3_2_minus[(j-3)*3]+secondOrderEdge3_3_minus[(j-3)*3];
                            secondOrderEdge3_origin_minus[(j-3)*3+1] = 2*secondOrderEdge3_minus[(j-3)*3+1]-3*secondOrderEdge3_2_minus[(j-3)*3+1]+secondOrderEdge3_3_minus[(j-3)*3+1];
                            secondOrderEdge3_origin_minus[(j-3)*3+2] = 2*secondOrderEdge3_minus[(j-3)*3+2]-3*secondOrderEdge3_2_minus[(j-3)*3+2]+secondOrderEdge3_3_minus[(j-3)*3+2];
                        }
                    }
                }
            }
            // rigid transformation
            float* firstOrderEdge4_minus_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge4_second_minus_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge4_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge4_second_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge2_minus_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge2_second_minus_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge2_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge2_second_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge1_minus_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge1_second_minus_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge1_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge1_second_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge3_minus_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge3_second_minus_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge3_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* firstOrderEdge3_second_temp = (float*)malloc((numberControlPointUDirection-4)*3*sizeof(float));
            float* secondOrderEdge1_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge1_2_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge1_3_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge1_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge1_2_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge1_3_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge4_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge4_2_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge4_3_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge4_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge4_2_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge4_3_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge2_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge2_2_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge2_3_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge2_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge2_2_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge2_3_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge3_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge3_2_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge3_3_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge3_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge3_2_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            float* secondOrderEdge3_3_minus_temp = (float*)malloc((numberControlPointUDirection-6)*3*sizeof(float));
            Multiply2(firstOrderEdge4_temp, firstOrderEdge4, T1d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge4_second_temp, firstOrderEdge4_second, T1d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge4_minus_temp, firstOrderEdge4_minus, T1d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge4_second_minus_temp, firstOrderEdge4_second_minus, T1d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge2_temp, firstOrderEdge2, T2d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge2_second_temp, firstOrderEdge2_second, T2d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge2_minus_temp, firstOrderEdge2_minus, T2d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge2_second_minus_temp, firstOrderEdge2_second_minus, T2d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge3_temp, firstOrderEdge3, T3d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge3_second_temp, firstOrderEdge3_second, T3d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge3_minus_temp, firstOrderEdge3_minus, T3d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge3_second_minus_temp, firstOrderEdge3_second_minus, T3d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge1_temp, firstOrderEdge1, T4d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge1_second_temp, firstOrderEdge1_second, T4d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge1_minus_temp, firstOrderEdge1_minus, T4d, numberControlPointUDirection-4);
            Multiply2(firstOrderEdge1_second_minus_temp, firstOrderEdge1_second_minus, T4d, numberControlPointUDirection-4);
            Multiply2(secondOrderEdge1_temp, secondOrderEdge1, T4d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge1_2_temp, secondOrderEdge1_2, T4d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge1_3_temp, secondOrderEdge1_3, T4d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge1_minus_temp, secondOrderEdge1_minus, T4d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge1_2_minus_temp, secondOrderEdge1_2_minus, T4d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge1_3_minus_temp, secondOrderEdge1_3_minus, T4d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge4_temp, secondOrderEdge4, T1d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge4_2_temp, secondOrderEdge4_2, T1d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge4_3_temp, secondOrderEdge4_3, T1d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge4_minus_temp, secondOrderEdge4_minus, T1d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge4_2_minus_temp, secondOrderEdge4_2_minus, T1d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge4_3_minus_temp, secondOrderEdge4_3_minus, T1d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge2_temp, secondOrderEdge2, T2d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge2_2_temp, secondOrderEdge2_2, T2d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge2_3_temp, secondOrderEdge2_3, T2d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge2_minus_temp, secondOrderEdge2_minus, T2d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge2_2_minus_temp, secondOrderEdge2_2_minus, T2d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge2_3_minus_temp, secondOrderEdge2_3_minus, T2d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge3_temp, secondOrderEdge3, T3d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge3_2_temp, secondOrderEdge3_2, T3d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge3_3_temp, secondOrderEdge3_3, T3d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge3_minus_temp, secondOrderEdge3_minus, T3d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge3_2_minus_temp, secondOrderEdge3_2_minus, T3d, numberControlPointUDirection-6);
            Multiply2(secondOrderEdge3_3_minus_temp, secondOrderEdge3_3_minus, T3d, numberControlPointUDirection-6);

            // calculate the vectors after transformation
            for (int i=0;i<numberControlPointUDirection-4;i++) {
                firstOrderEdge1_second_origin[i*3] = (-1)*(firstOrderEdge1_temp[i*3]-firstOrderEdge1_second_temp[i*3]);
                firstOrderEdge1_second_origin[i*3+1] = (-1)*(firstOrderEdge1_temp[i*3+1]-firstOrderEdge1_second_temp[i*3+1]);
                firstOrderEdge1_second_origin[i*3+2] = (-1)*(firstOrderEdge1_temp[i*3+2]-firstOrderEdge1_second_temp[i*3+2]);
                firstOrderEdge1_second_origin_minus[i*3] = (-1)*(firstOrderEdge1_minus_temp[i*3]-firstOrderEdge1_second_minus_temp[i*3]);
                firstOrderEdge1_second_origin_minus[i*3+1] = (-1)*(firstOrderEdge1_minus_temp[i*3+1]-firstOrderEdge1_second_minus_temp[i*3+1]);
                firstOrderEdge1_second_origin_minus[i*3+2] = (-1)*(firstOrderEdge1_minus_temp[i*3+2]-firstOrderEdge1_second_minus_temp[i*3+2]);
                firstOrderEdge2_second_origin[i*3] = (-1)*(firstOrderEdge2_temp[i*3]-firstOrderEdge2_second_temp[i*3]);
                firstOrderEdge2_second_origin[i*3+1] = (-1)*(firstOrderEdge2_temp[i*3+1]-firstOrderEdge2_second_temp[i*3+1]);
                firstOrderEdge2_second_origin[i*3+2] = (-1)*(firstOrderEdge2_temp[i*3+2]-firstOrderEdge2_second_temp[i*3+2]);
                firstOrderEdge2_second_origin_minus[i*3] = (-1)*(firstOrderEdge2_minus_temp[i*3]-firstOrderEdge2_second_minus_temp[i*3]);
                firstOrderEdge2_second_origin_minus[i*3+1] = (-1)*(firstOrderEdge2_minus_temp[i*3+1]-firstOrderEdge2_second_minus_temp[i*3+1]);
                firstOrderEdge2_second_origin_minus[i*3+2] = (-1)*(firstOrderEdge2_minus_temp[i*3+2]-firstOrderEdge2_second_minus_temp[i*3+2]);
                firstOrderEdge3_second_origin[i*3] = (-1)*(firstOrderEdge3_temp[i*3]-firstOrderEdge3_second_temp[i*3]);
                firstOrderEdge3_second_origin[i*3+1] = (-1)*(firstOrderEdge3_temp[i*3+1]-firstOrderEdge3_second_temp[i*3+1]);
                firstOrderEdge3_second_origin[i*3+2] = (-1)*(firstOrderEdge3_temp[i*3+2]-firstOrderEdge3_second_temp[i*3+2]);
                firstOrderEdge3_second_origin_minus[i*3] = (-1)*(firstOrderEdge3_minus_temp[i*3]-firstOrderEdge3_second_minus_temp[i*3]);
                firstOrderEdge3_second_origin_minus[i*3+1] = (-1)*(firstOrderEdge3_minus_temp[i*3+1]-firstOrderEdge3_second_minus_temp[i*3+1]);
                firstOrderEdge3_second_origin_minus[i*3+2] = (-1)*(firstOrderEdge3_minus_temp[i*3+2]-firstOrderEdge3_second_minus_temp[i*3+2]);
                firstOrderEdge4_second_origin[i*3] = (-1)*(firstOrderEdge4_temp[i*3]-firstOrderEdge4_second_temp[i*3]);
                firstOrderEdge4_second_origin[i*3+1] = (-1)*(firstOrderEdge4_temp[i*3+1]-firstOrderEdge4_second_temp[i*3+1]);
                firstOrderEdge4_second_origin[i*3+2] = (-1)*(firstOrderEdge4_temp[i*3+2]-firstOrderEdge4_second_temp[i*3+2]);
                firstOrderEdge4_second_origin_minus[i*3] = (-1)*(firstOrderEdge4_minus_temp[i*3]-firstOrderEdge4_second_minus_temp[i*3]);
                firstOrderEdge4_second_origin_minus[i*3+1] = (-1)*(firstOrderEdge4_minus_temp[i*3+1]-firstOrderEdge4_second_minus_temp[i*3+1]);
                firstOrderEdge4_second_origin_minus[i*3+2] = (-1)*(firstOrderEdge4_minus_temp[i*3+2]-firstOrderEdge4_second_minus_temp[i*3+2]);
            }
            for (int i=0;i<numberControlPointUDirection-6;i++) {
                secondOrderEdge1_second_origin[i*3] = (-1)*(2*secondOrderEdge1_temp[i*3]-3*secondOrderEdge1_2_temp[i*3]+secondOrderEdge1_3_temp[i*3]);
                secondOrderEdge1_second_origin[i*3+1] = (-1)*(2*secondOrderEdge1_temp[i*3+1]-3*secondOrderEdge1_2_temp[i*3+1]+secondOrderEdge1_3_temp[i*3+1]);
                secondOrderEdge1_second_origin[i*3+2] = (-1)*(2*secondOrderEdge1_temp[i*3+2]-3*secondOrderEdge1_2_temp[i*3+2]+secondOrderEdge1_3_temp[i*3+2]);
                secondOrderEdge2_second_origin[i*3] = (-1)*(2*secondOrderEdge2_temp[i*3]-3*secondOrderEdge2_2_temp[i*3]+secondOrderEdge2_3_temp[i*3]);
                secondOrderEdge2_second_origin[i*3+1] = (-1)*(2*secondOrderEdge2_temp[i*3+1]-3*secondOrderEdge2_2_temp[i*3+1]+secondOrderEdge2_3_temp[i*3+1]);
                secondOrderEdge2_second_origin[i*3+2] = (-1)*(2*secondOrderEdge2_temp[i*3+2]-3*secondOrderEdge2_2_temp[i*3+2]+secondOrderEdge2_3_temp[i*3+2]);
                secondOrderEdge3_second_origin[i*3] = (-1)*(2*secondOrderEdge3_temp[i*3]-3*secondOrderEdge3_2_temp[i*3]+secondOrderEdge3_3_temp[i*3]);
                secondOrderEdge3_second_origin[i*3+1] = (-1)*(2*secondOrderEdge3_temp[i*3+1]-3*secondOrderEdge3_2_temp[i*3+1]+secondOrderEdge3_3_temp[i*3+1]);
                secondOrderEdge3_second_origin[i*3+2] = (-1)*(2*secondOrderEdge3_temp[i*3+2]-3*secondOrderEdge3_2_temp[i*3+2]+secondOrderEdge3_3_temp[i*3+2]);
                secondOrderEdge4_second_origin[i*3] = (-1)*(2*secondOrderEdge4_temp[i*3]-3*secondOrderEdge4_2_temp[i*3]+secondOrderEdge4_3_temp[i*3]);
                secondOrderEdge4_second_origin[i*3+1] = (-1)*(2*secondOrderEdge4_temp[i*3+1]-3*secondOrderEdge4_2_temp[i*3+1]+secondOrderEdge4_3_temp[i*3+1]);
                secondOrderEdge4_second_origin[i*3+2] = (-1)*(2*secondOrderEdge4_temp[i*3+2]-3*secondOrderEdge4_2_temp[i*3+2]+secondOrderEdge4_3_temp[i*3+2]);
                secondOrderEdge1_second_origin_minus[i*3] = (-1)*(2*secondOrderEdge1_minus_temp[i*3]-3*secondOrderEdge1_2_minus_temp[i*3]+secondOrderEdge1_3_minus_temp[i*3]);
                secondOrderEdge1_second_origin_minus[i*3+1] = (-1)*(2*secondOrderEdge1_minus_temp[i*3+1]-3*secondOrderEdge1_2_minus_temp[i*3+1]+secondOrderEdge1_3_minus_temp[i*3+1]);
                secondOrderEdge1_second_origin_minus[i*3+2] = (-1)*(2*secondOrderEdge1_minus_temp[i*3+2]-3*secondOrderEdge1_2_minus_temp[i*3+2]+secondOrderEdge1_3_minus_temp[i*3+2]);
                secondOrderEdge2_second_origin_minus[i*3] = (-1)*(2*secondOrderEdge2_minus_temp[i*3]-3*secondOrderEdge2_2_minus_temp[i*3]+secondOrderEdge2_3_minus_temp[i*3]);
                secondOrderEdge2_second_origin_minus[i*3+1] = (-1)*(2*secondOrderEdge2_minus_temp[i*3+1]-3*secondOrderEdge2_2_minus_temp[i*3+1]+secondOrderEdge2_3_minus_temp[i*3+1]);
                secondOrderEdge2_second_origin_minus[i*3+2] = (-1)*(2*secondOrderEdge2_minus_temp[i*3+2]-3*secondOrderEdge2_2_minus_temp[i*3+2]+secondOrderEdge2_3_minus_temp[i*3+2]);
                secondOrderEdge3_second_origin_minus[i*3] = (-1)*(2*secondOrderEdge3_minus_temp[i*3]-3*secondOrderEdge3_2_minus_temp[i*3]+secondOrderEdge3_3_minus_temp[i*3]);
                secondOrderEdge3_second_origin_minus[i*3+1] = (-1)*(2*secondOrderEdge3_minus_temp[i*3+1]-3*secondOrderEdge3_2_minus_temp[i*3+1]+secondOrderEdge3_3_minus_temp[i*3+1]);
                secondOrderEdge3_second_origin_minus[i*3+2] = (-1)*(2*secondOrderEdge3_minus_temp[i*3+2]-3*secondOrderEdge3_2_minus_temp[i*3+2]+secondOrderEdge3_3_minus_temp[i*3+2]);
                secondOrderEdge4_second_origin_minus[i*3] = (-1)*(2*secondOrderEdge4_minus_temp[i*3]-3*secondOrderEdge4_2_minus_temp[i*3]+secondOrderEdge4_3_minus_temp[i*3]);
                secondOrderEdge4_second_origin_minus[i*3+1] = (-1)*(2*secondOrderEdge4_minus_temp[i*3+1]-3*secondOrderEdge4_2_minus_temp[i*3+1]+secondOrderEdge4_3_minus_temp[i*3+1]);
                secondOrderEdge4_second_origin_minus[i*3+2] = (-1)*(2*secondOrderEdge4_minus_temp[i*3+2]-3*secondOrderEdge4_2_minus_temp[i*3+2]+secondOrderEdge4_3_minus_temp[i*3+2]);
            }

            // calculate the constraints for the plus surface
            for (int i=0;i<numberControlPointUDirection;i++) {
                for (int j=0;j<numberControlPointVDirection;j++) {
                    // edge 1
                    if (i == numberControlPointUDirection-2 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge4_second_origin_minus[(j-2)*3]-firstOrderEdge1_origin[(j-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge4_second_origin_minus[(j-2)*3+1]-firstOrderEdge1_origin[(j-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge4_second_origin_minus[(j-2)*3+2]-firstOrderEdge1_origin[(j-2)*3+2])/2;
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3] = (secondOrderEdge4_second_origin_minus[(j-3)*3]-secondOrderEdge1_origin[(j-3)*3])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3+1] = (secondOrderEdge4_second_origin_minus[(j-3)*3+1]-secondOrderEdge1_origin[(j-3)*3+1])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3+2] = (secondOrderEdge4_second_origin_minus[(j-3)*3+2]-secondOrderEdge1_origin[(j-3)*3+2])/12.0;
                        }
                    }
                    // edge 2
                    if (j == 1 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge2_second_origin_minus[(i-2)*3]-firstOrderEdge2_origin[(i-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge2_second_origin_minus[(i-2)*3+1]-firstOrderEdge2_origin[(i-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge2_second_origin_minus[(i-2)*3+2]-firstOrderEdge2_origin[(i-2)*3+2])/2;
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3] = (secondOrderEdge2_second_origin_minus[(i-3)*3]-secondOrderEdge2_origin[(i-3)*3])/12.0;
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3+1] = (secondOrderEdge2_second_origin_minus[(i-3)*3+1]-secondOrderEdge2_origin[(i-3)*3+1])/12.0;
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3+2] = (secondOrderEdge2_second_origin_minus[(i-3)*3+2]-secondOrderEdge2_origin[(i-3)*3+2])/12.0;
                        }
                    }
                    // edge 3
                    if (i == 1 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge3_second_origin_minus[(j-2)*3]-firstOrderEdge3_origin[(j-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge3_second_origin_minus[(j-2)*3+1]-firstOrderEdge3_origin[(j-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge3_second_origin_minus[(j-2)*3+2]-firstOrderEdge3_origin[(j-2)*3+2])/2;
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3] = (secondOrderEdge3_second_origin_minus[(j-3)*3]-secondOrderEdge3_origin[(j-3)*3])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3+1] = (secondOrderEdge3_second_origin_minus[(j-3)*3+1]-secondOrderEdge3_origin[(j-3)*3+1])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3+2] = (secondOrderEdge3_second_origin_minus[(j-3)*3+2]-secondOrderEdge3_origin[(j-3)*3+2])/12.0;
                        }
                    }
                    // edge 4
                    if (j == numberControlPointVDirection-2 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge1_second_origin_minus[(i-2)*3]-firstOrderEdge4_origin[(i-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge1_second_origin_minus[(i-2)*3+1]-firstOrderEdge4_origin[(i-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge1_second_origin_minus[(i-2)*3+2]-firstOrderEdge4_origin[(i-2)*3+2])/2;
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3] = (secondOrderEdge1_second_origin_minus[(i-3)*3]-secondOrderEdge4_origin[(i-3)*3])/12.0;
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3+1] = (secondOrderEdge1_second_origin_minus[(i-3)*3+1]-secondOrderEdge4_origin[(i-3)*3+1])/12.0;
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3+2] = (secondOrderEdge1_second_origin_minus[(i-3)*3+2]-secondOrderEdge4_origin[(i-3)*3+2])/12.0;
                        }
                    }
                }
            }
        }
        else if (modelType == "SchwarzP") {
            for (int i=0;i<numberControlPointUDirection;i++) {
                for (int j=0;j<numberControlPointVDirection;j++) {
                    // edge 1
                    if (i == 0 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        firstOrderEdge1[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge1[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge1[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge1_second[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3];
                        firstOrderEdge1_second[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+1];
                        firstOrderEdge1_second[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+2];
                        firstOrderEdge1_origin[(j-2)*3] = firstOrderEdge1[(j-2)*3]-firstOrderEdge1_second[(j-2)*3];
                        firstOrderEdge1_origin[(j-2)*3+1] = firstOrderEdge1[(j-2)*3+1]-firstOrderEdge1_second[(j-2)*3+1];
                        firstOrderEdge1_origin[(j-2)*3+2] = firstOrderEdge1[(j-2)*3+2]-firstOrderEdge1_second[(j-2)*3+2];
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            secondOrderEdge1[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge1[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge1[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge1_2[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3];
                            secondOrderEdge1_2[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+1];
                            secondOrderEdge1_2[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+1))*3+2];
                            secondOrderEdge1_3[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3];
                            secondOrderEdge1_3[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3+1];
                            secondOrderEdge1_3[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i+2))*3+2];
                            secondOrderEdge1_origin[(j-3)*3] = 2*secondOrderEdge1[(j-3)*3]-3*secondOrderEdge1_2[(j-3)*3]+secondOrderEdge1_3[(j-3)*3];
                            secondOrderEdge1_origin[(j-3)*3+1] = 2*secondOrderEdge1[(j-3)*3+1]-3*secondOrderEdge1_2[(j-3)*3+1]+secondOrderEdge1_3[(j-3)*3+1];
                            secondOrderEdge1_origin[(j-3)*3+2] = 2*secondOrderEdge1[(j-3)*3+2]-3*secondOrderEdge1_2[(j-3)*3+2]+secondOrderEdge1_3[(j-3)*3+2];
                        }
                    }
                    // edge 2
                    if (j == 0 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        firstOrderEdge2[(i-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge2[(i-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge2[(i-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge2_second[(i-2)*3] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3];
                        firstOrderEdge2_second[(i-2)*3+1] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge2_second[(i-2)*3+2] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge2_origin[(i-2)*3] = firstOrderEdge2[(i-2)*3]-firstOrderEdge2_second[(i-2)*3];
                        firstOrderEdge2_origin[(i-2)*3+1] = firstOrderEdge2[(i-2)*3+1]-firstOrderEdge2_second[(i-2)*3+1];
                        firstOrderEdge2_origin[(i-2)*3+2] = firstOrderEdge2[(i-2)*3+2]-firstOrderEdge2_second[(i-2)*3+2];
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            secondOrderEdge2[(i-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge2[(i-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2[(i-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_2[(i-3)*3] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3];
                            secondOrderEdge2_2[(i-3)*3+1] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2_2[(i-3)*3+2] = controlPointsPlusPIAGPU[((j+1)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_3[(i-3)*3] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3];
                            secondOrderEdge2_3[(i-3)*3+1] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge2_3[(i-3)*3+2] = controlPointsPlusPIAGPU[((j+2)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge2_origin[(i-3)*3] = 2*secondOrderEdge2[(i-3)*3]-3*secondOrderEdge2_2[(i-3)*3]+secondOrderEdge2_3[(i-3)*3];
                            secondOrderEdge2_origin[(i-3)*3+1] = 2*secondOrderEdge2[(i-3)*3+1]-3*secondOrderEdge2_2[(i-3)*3+1]+secondOrderEdge2_3[(i-3)*3+1];
                            secondOrderEdge2_origin[(i-3)*3+2] = 2*secondOrderEdge2[(i-3)*3+2]-3*secondOrderEdge2_2[(i-3)*3+2]+secondOrderEdge2_3[(i-3)*3+2];
                        }
                    }
                    // edge 3
                    if (i == numberControlPointUDirection-1 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        firstOrderEdge3[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge3[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge3[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge3_second[(j-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3];
                        firstOrderEdge3_second[(j-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+1];
                        firstOrderEdge3_second[(j-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+2];
                        firstOrderEdge3_origin[(j-2)*3] = firstOrderEdge3[(j-2)*3]-firstOrderEdge3_second[(j-2)*3];
                        firstOrderEdge3_origin[(j-2)*3+1] = firstOrderEdge3[(j-2)*3+1]-firstOrderEdge3_second[(j-2)*3+1];
                        firstOrderEdge3_origin[(j-2)*3+2] = firstOrderEdge3[(j-2)*3+2]-firstOrderEdge3_second[(j-2)*3+2];
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            secondOrderEdge3[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge3[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge3[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge3_2[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3];
                            secondOrderEdge3_2[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+1];
                            secondOrderEdge3_2[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-1))*3+2];
                            secondOrderEdge3_3[(j-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3];
                            secondOrderEdge3_3[(j-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3+1];
                            secondOrderEdge3_3[(j-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+(i-2))*3+2];
                            secondOrderEdge3_origin[(j-3)*3] = 2*secondOrderEdge3[(j-3)*3]-3*secondOrderEdge3_2[(j-3)*3]+secondOrderEdge3_3[(j-3)*3];
                            secondOrderEdge3_origin[(j-3)*3+1] = 2*secondOrderEdge3[(j-3)*3+1]-3*secondOrderEdge3_2[(j-3)*3+1]+secondOrderEdge3_3[(j-3)*3+1];
                            secondOrderEdge3_origin[(j-3)*3+2] = 2*secondOrderEdge3[(j-3)*3+2]-3*secondOrderEdge3_2[(j-3)*3+2]+secondOrderEdge3_3[(j-3)*3+2];
                        }
                    }
                    // edge 4
                    if (j == numberControlPointVDirection-1 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        firstOrderEdge4[(i-2)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                        firstOrderEdge4[(i-2)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge4[(i-2)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge4_second[(i-2)*3] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3];
                        firstOrderEdge4_second[(i-2)*3+1] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+1];
                        firstOrderEdge4_second[(i-2)*3+2] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+2];
                        firstOrderEdge4_origin[(i-2)*3] = firstOrderEdge4[(i-2)*3]-firstOrderEdge4_second[(i-2)*3];
                        firstOrderEdge4_origin[(i-2)*3+1] = firstOrderEdge4[(i-2)*3+1]-firstOrderEdge4_second[(i-2)*3+1];
                        firstOrderEdge4_origin[(i-2)*3+2] = firstOrderEdge4[(i-2)*3+2]-firstOrderEdge4_second[(i-2)*3+2];
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            secondOrderEdge4[(i-3)*3] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3];
                            secondOrderEdge4[(i-3)*3+1] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4[(i-3)*3+2] = controlPointsPlusPIAGPU[(j*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_2[(i-3)*3] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3];
                            secondOrderEdge4_2[(i-3)*3+1] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4_2[(i-3)*3+2] = controlPointsPlusPIAGPU[((j-1)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_3[(i-3)*3] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3];
                            secondOrderEdge4_3[(i-3)*3+1] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3+1];
                            secondOrderEdge4_3[(i-3)*3+2] = controlPointsPlusPIAGPU[((j-2)*numberControlPointUDirection+i)*3+2];
                            secondOrderEdge4_origin[(i-3)*3] = 2*secondOrderEdge4[(i-3)*3]-3*secondOrderEdge4_2[(i-3)*3]+secondOrderEdge4_3[(i-3)*3];
                            secondOrderEdge4_origin[(i-3)*3+1] = 2*secondOrderEdge4[(i-3)*3+1]-3*secondOrderEdge4_2[(i-3)*3+1]+secondOrderEdge4_3[(i-3)*3+1];
                            secondOrderEdge4_origin[(i-3)*3+2] = 2*secondOrderEdge4[(i-3)*3+2]-3*secondOrderEdge4_2[(i-3)*3+2]+secondOrderEdge4_3[(i-3)*3+2];
                        }
                    }
                }
            }
            // rigid transformation
            Multiply(firstOrderEdge1, T1p, numberControlPointUDirection-4);
            Multiply(firstOrderEdge1_second, T1p, numberControlPointUDirection-4);
            Multiply(secondOrderEdge1, T1p, numberControlPointUDirection-6);
            Multiply(secondOrderEdge1_2, T1p, numberControlPointUDirection-6);
            Multiply(secondOrderEdge1_3, T1p, numberControlPointUDirection-6);
            Multiply(firstOrderEdge2, T2p, numberControlPointUDirection-4);
            Multiply(firstOrderEdge2_second, T2p, numberControlPointUDirection-4);
            Multiply(secondOrderEdge2, T2p, numberControlPointUDirection-6);
            Multiply(secondOrderEdge2_2, T2p, numberControlPointUDirection-6);
            Multiply(secondOrderEdge2_3, T2p, numberControlPointUDirection-6);
            Multiply(firstOrderEdge3, T3p, numberControlPointUDirection-4);
            Multiply(firstOrderEdge3_second, T3p, numberControlPointUDirection-4);
            Multiply(secondOrderEdge3, T3p, numberControlPointUDirection-6);
            Multiply(secondOrderEdge3_2, T3p, numberControlPointUDirection-6);
            Multiply(secondOrderEdge3_3, T3p, numberControlPointUDirection-6);
            Multiply(firstOrderEdge4, T4p, numberControlPointUDirection-4);
            Multiply(firstOrderEdge4_second, T4p, numberControlPointUDirection-4);
            Multiply(secondOrderEdge4, T4p, numberControlPointUDirection-6);
            Multiply(secondOrderEdge4_2, T4p, numberControlPointUDirection-6);
            Multiply(secondOrderEdge4_3, T4p, numberControlPointUDirection-6);

            for (int i=0;i<numberControlPointUDirection-4;i++) {
                firstOrderEdge1_second_origin[i*3] = (-1)*(firstOrderEdge1[i*3]-firstOrderEdge1_second[i*3]);
                firstOrderEdge1_second_origin[i*3+1] = (-1)*(firstOrderEdge1[i*3+1]-firstOrderEdge1_second[i*3+1]);
                firstOrderEdge1_second_origin[i*3+2] = (-1)*(firstOrderEdge1[i*3+2]-firstOrderEdge1_second[i*3+2]);
                firstOrderEdge2_second_origin[i*3] = (-1)*(firstOrderEdge2[i*3]-firstOrderEdge2_second[i*3]);
                firstOrderEdge2_second_origin[i*3+1] = (-1)*(firstOrderEdge2[i*3+1]-firstOrderEdge2_second[i*3+1]);
                firstOrderEdge2_second_origin[i*3+2] = (-1)*(firstOrderEdge2[i*3+2]-firstOrderEdge2_second[i*3+2]);
                firstOrderEdge3_second_origin[i*3] = (-1)*(firstOrderEdge3[i*3]-firstOrderEdge3_second[i*3]);
                firstOrderEdge3_second_origin[i*3+1] = (-1)*(firstOrderEdge3[i*3+1]-firstOrderEdge3_second[i*3+1]);
                firstOrderEdge3_second_origin[i*3+2] = (-1)*(firstOrderEdge3[i*3+2]-firstOrderEdge3_second[i*3+2]);
                firstOrderEdge4_second_origin[i*3] = (-1)*(firstOrderEdge4[i*3]-firstOrderEdge4_second[i*3]);
                firstOrderEdge4_second_origin[i*3+1] = (-1)*(firstOrderEdge4[i*3+1]-firstOrderEdge4_second[i*3+1]);
                firstOrderEdge4_second_origin[i*3+2] = (-1)*(firstOrderEdge4[i*3+2]-firstOrderEdge4_second[i*3+2]);
            }
            for (int i=0;i<numberControlPointUDirection-6;i++) {
                secondOrderEdge1_second_origin[i*3] = (-1)*(2*secondOrderEdge1[i*3]-3*secondOrderEdge1_2[i*3]+secondOrderEdge1_3[i*3]);
                secondOrderEdge1_second_origin[i*3+1] = (-1)*(2*secondOrderEdge1[i*3+1]-3*secondOrderEdge1_2[i*3+1]+secondOrderEdge1_3[i*3+1]);
                secondOrderEdge1_second_origin[i*3+2] = (-1)*(2*secondOrderEdge1[i*3+2]-3*secondOrderEdge1_2[i*3+2]+secondOrderEdge1_3[i*3+2]);
                secondOrderEdge2_second_origin[i*3] = (-1)*(2*secondOrderEdge2[i*3]-3*secondOrderEdge2_2[i*3]+secondOrderEdge2_3[i*3]);
                secondOrderEdge2_second_origin[i*3+1] = (-1)*(2*secondOrderEdge2[i*3+1]-3*secondOrderEdge2_2[i*3+1]+secondOrderEdge2_3[i*3+1]);
                secondOrderEdge2_second_origin[i*3+2] = (-1)*(2*secondOrderEdge2[i*3+2]-3*secondOrderEdge2_2[i*3+2]+secondOrderEdge2_3[i*3+2]);
                secondOrderEdge3_second_origin[i*3] = (-1)*(2*secondOrderEdge3[i*3]-3*secondOrderEdge3_2[i*3]+secondOrderEdge3_3[i*3]);
                secondOrderEdge3_second_origin[i*3+1] = (-1)*(2*secondOrderEdge3[i*3+1]-3*secondOrderEdge3_2[i*3+1]+secondOrderEdge3_3[i*3+1]);
                secondOrderEdge3_second_origin[i*3+2] = (-1)*(2*secondOrderEdge3[i*3+2]-3*secondOrderEdge3_2[i*3+2]+secondOrderEdge3_3[i*3+2]);
                secondOrderEdge4_second_origin[i*3] = (-1)*(2*secondOrderEdge4[i*3]-3*secondOrderEdge4_2[i*3]+secondOrderEdge4_3[i*3]);
                secondOrderEdge4_second_origin[i*3+1] = (-1)*(2*secondOrderEdge4[i*3+1]-3*secondOrderEdge4_2[i*3+1]+secondOrderEdge4_3[i*3+1]);
                secondOrderEdge4_second_origin[i*3+2] = (-1)*(2*secondOrderEdge4[i*3+2]-3*secondOrderEdge4_2[i*3+2]+secondOrderEdge4_3[i*3+2]);
            }

            for (int i=0;i<numberControlPointUDirection;i++) {
                for (int j=0;j<numberControlPointVDirection;j++) {
                    // edge 1
                    if (i == 1 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge1_second_origin[(numberControlPointUDirection-5-(j-2))*3]-firstOrderEdge1_origin[(j-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge1_second_origin[(numberControlPointUDirection-5-(j-2))*3+1]-firstOrderEdge1_origin[(j-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge1_second_origin[(numberControlPointUDirection-5-(j-2))*3+2]-firstOrderEdge1_origin[(j-2)*3+2])/2;
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3] = (secondOrderEdge1_second_origin[(j-3)*3]-secondOrderEdge1_origin[(j-3)*3])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3+1] = (secondOrderEdge1_second_origin[(j-3)*3+1]-secondOrderEdge1_origin[(j-3)*3+1])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3+2] = (secondOrderEdge1_second_origin[(j-3)*3+2]-secondOrderEdge1_origin[(j-3)*3+2])/12.0;
                        }
                    }
                    // edge 2
                    if (j == 1 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge2_second_origin[(i-2)*3]-firstOrderEdge2_origin[(i-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge2_second_origin[(i-2)*3+1]-firstOrderEdge2_origin[(i-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge2_second_origin[(i-2)*3+2]-firstOrderEdge2_origin[(i-2)*3+2])/2;
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3] = (secondOrderEdge2_second_origin[(i-3)*3]-secondOrderEdge2_origin[(i-3)*3])/12.0;
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3+1] = (secondOrderEdge2_second_origin[(i-3)*3+1]-secondOrderEdge2_origin[(i-3)*3+1])/12.0;
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3+2] = (secondOrderEdge2_second_origin[(i-3)*3+2]-secondOrderEdge2_origin[(i-3)*3+2])/12.0;
                        }
                    }
                    // edge 3
                    if (i == numberControlPointUDirection-2 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge3_second_origin[(j-2)*3]-firstOrderEdge3_origin[(j-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge3_second_origin[(j-2)*3+1]-firstOrderEdge3_origin[(j-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge3_second_origin[(j-2)*3+2]-firstOrderEdge3_origin[(j-2)*3+2])/2;
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3] = (secondOrderEdge3_second_origin[(j-3)*3]-secondOrderEdge3_origin[(j-3)*3])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3+1] = (secondOrderEdge3_second_origin[(j-3)*3+1]-secondOrderEdge3_origin[(j-3)*3+1])/12.0;
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3+2] = (secondOrderEdge3_second_origin[(j-3)*3+2]-secondOrderEdge3_origin[(j-3)*3+2])/12.0;
                        }
                    }
                    // edge 4
                    if (j == numberControlPointVDirection-2 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge4_second_origin[(i-2)*3]-firstOrderEdge4_origin[(i-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge4_second_origin[(i-2)*3+1]-firstOrderEdge4_origin[(i-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge4_second_origin[(i-2)*3+2]-firstOrderEdge4_origin[(i-2)*3+2])/2;
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3] = (secondOrderEdge4_second_origin[(i-3)*3]-secondOrderEdge4_origin[(i-3)*3])/12.0;
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3+1] = (secondOrderEdge4_second_origin[(i-3)*3+1]-secondOrderEdge4_origin[(i-3)*3+1])/12.0;
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3+2] = (secondOrderEdge4_second_origin[(i-3)*3+2]-secondOrderEdge4_origin[(i-3)*3+2])/12.0;
                        }
                    }
                }
            }
        }
        
        if (count>=constraintIterationNumber) {
            for (int i=0;i<numberControlPointUDirection*numberControlPointVDirection;i++) {
                controlPointsPlusPIAGPU[i*3] = controlPointsPlusPIAGPU[i*3] + surfacePlusMatrix(i,0) - evaluationSurfacePlusPIA(i,0) + constrainMatrix[i*3];
                controlPointsPlusPIAGPU[i*3+1] = controlPointsPlusPIAGPU[i*3+1] + surfacePlusMatrix(i,1) - evaluationSurfacePlusPIA(i,1) + constrainMatrix[i*3+1];
                controlPointsPlusPIAGPU[i*3+2] = controlPointsPlusPIAGPU[i*3+2] + surfacePlusMatrix(i,2) - evaluationSurfacePlusPIA(i,2) + constrainMatrix[i*3+2];
            }
        } else {
            for (int i=0;i<numberControlPointUDirection*numberControlPointVDirection;i++) {
                controlPointsPlusPIAGPU[i*3] = controlPointsPlusPIAGPU[i*3] + surfacePlusMatrix(i,0) - evaluationSurfacePlusPIA(i,0);
                controlPointsPlusPIAGPU[i*3+1] = controlPointsPlusPIAGPU[i*3+1] + surfacePlusMatrix(i,1) - evaluationSurfacePlusPIA(i,1);
                controlPointsPlusPIAGPU[i*3+2] = controlPointsPlusPIAGPU[i*3+2] + surfacePlusMatrix(i,2) - evaluationSurfacePlusPIA(i,2);
            }
        }
        if (modelType == "Gyroid") {
            clearMatrix(constrainMatrix, numberControlPointUDirection*numberControlPointVDirection);
        } else if (modelType == "Diamond") {
            clearMatrix(constrainMatrix, numberControlPointUDirection*numberControlPointVDirection);
            // calculate the constraints for minus surface
            for (int i=0;i<numberControlPointUDirection;i++) {
                for (int j=0;j<numberControlPointVDirection;j++) {
                    // edge 1
                    if (i == numberControlPointUDirection-2 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge4_second_origin[(j-2)*3]-firstOrderEdge1_origin_minus[(j-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge4_second_origin[(j-2)*3+1]-firstOrderEdge1_origin_minus[(j-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge4_second_origin[(j-2)*3+2]-firstOrderEdge1_origin_minus[(j-2)*3+2])/2;
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3] = (secondOrderEdge4_second_origin[(j-3)*3]-secondOrderEdge1_origin_minus[(j-3)*3])/2;
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3+1] = (secondOrderEdge4_second_origin[(j-3)*3+1]-secondOrderEdge1_origin_minus[(j-3)*3+1])/2;
                            constrainMatrix[(j*numberControlPointUDirection+(i-1))*3+2] = (secondOrderEdge4_second_origin[(j-3)*3+2]-secondOrderEdge1_origin_minus[(j-3)*3+2])/2;
                        }
                    }
                    // edge 2
                    if (j == 1 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge2_second_origin[(i-2)*3]-firstOrderEdge2_origin_minus[(i-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge2_second_origin[(i-2)*3+1]-firstOrderEdge2_origin_minus[(i-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge2_second_origin[(i-2)*3+2]-firstOrderEdge2_origin_minus[(i-2)*3+2])/2;
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3] = (secondOrderEdge2_second_origin[(i-3)*3]-secondOrderEdge2_origin_minus[(i-3)*3])/2;
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3+1] = (secondOrderEdge2_second_origin[(i-3)*3+1]-secondOrderEdge2_origin_minus[(i-3)*3+1])/2;
                            constrainMatrix[((j+1)*numberControlPointUDirection+i)*3+2] = (secondOrderEdge2_second_origin[(i-3)*3+2]-secondOrderEdge2_origin_minus[(i-3)*3+2])/2;
                        }
                    }
                    // edge 3
                    if (i == 1 && j != 0 && j != 1 && j != numberControlPointVDirection-1 && j != numberControlPointVDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge3_second_origin[(j-2)*3]-firstOrderEdge3_origin_minus[(j-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge3_second_origin[(j-2)*3+1]-firstOrderEdge3_origin_minus[(j-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge3_second_origin[(j-2)*3+2]-firstOrderEdge3_origin_minus[(j-2)*3+2])/2;
                        if (j != 2 && j != numberControlPointVDirection-3) {
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3] = (secondOrderEdge3_second_origin[(j-3)*3]-secondOrderEdge3_origin_minus[(j-3)*3])/2;
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3+1] = (secondOrderEdge3_second_origin[(j-3)*3+1]-secondOrderEdge3_origin_minus[(j-3)*3+1])/2;
                            constrainMatrix[(j*numberControlPointUDirection+(i+1))*3+2] = (secondOrderEdge3_second_origin[(j-3)*3+2]-secondOrderEdge3_origin_minus[(j-3)*3+2])/2;
                        }
                    }
                    // edge 4
                    if (j == numberControlPointVDirection-2 && i != 0 && i != 1 && i != numberControlPointUDirection-1 && i != numberControlPointUDirection-2) {
                        constrainMatrix[(j*numberControlPointUDirection+i)*3] = (-1)*(firstOrderEdge1_second_origin[(i-2)*3]-firstOrderEdge4_origin_minus[(i-2)*3])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+1] = (-1)*(firstOrderEdge1_second_origin[(i-2)*3+1]-firstOrderEdge4_origin_minus[(i-2)*3+1])/2;
                        constrainMatrix[(j*numberControlPointUDirection+i)*3+2] = (-1)*(firstOrderEdge1_second_origin[(i-2)*3+2]-firstOrderEdge4_origin_minus[(i-2)*3+2])/2;
                        if (i != 2 && i != numberControlPointUDirection-3) {
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3] = (secondOrderEdge1_second_origin[(i-3)*3]-secondOrderEdge4_origin_minus[(i-3)*3])/2;
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3+1] = (secondOrderEdge1_second_origin[(i-3)*3+1]-secondOrderEdge4_origin_minus[(i-3)*3+1])/2;
                            constrainMatrix[((j-1)*numberControlPointUDirection+i)*3+2] = (secondOrderEdge1_second_origin[(i-3)*3+2]-secondOrderEdge4_origin_minus[(i-3)*3+2])/2;
                        }
                    }
                }
            }
        }
        else if (modelType == "SchwarzP") {
            clearMatrix(constrainMatrix, numberControlPointUDirection*numberControlPointVDirection);
        }

        if (modelType != "SchwarzP") {
            if (count>=constraintIterationNumber) {
                for (int i=0;i<numberControlPointUDirection*numberControlPointVDirection;i++) {
                    controlPointsMinusPIAGPU[i*3] = controlPointsMinusPIAGPU[i*3] + surfaceMinusMatrix(i,0) - evaluationSurfaceMinusPIA(i,0) + constrainMatrix[i*3];
                    controlPointsMinusPIAGPU[i*3+1] = controlPointsMinusPIAGPU[i*3+1] + surfaceMinusMatrix(i,1) - evaluationSurfaceMinusPIA(i,1) + constrainMatrix[i*3+1];
                    controlPointsMinusPIAGPU[i*3+2] = controlPointsMinusPIAGPU[i*3+2] + surfaceMinusMatrix(i,2) - evaluationSurfaceMinusPIA(i,2) + constrainMatrix[i*3+2];
                }
            } else {
                for (int i=0;i<numberControlPointUDirection*numberControlPointVDirection;i++) {
                    controlPointsMinusPIAGPU[i*3] = controlPointsMinusPIAGPU[i*3] + surfaceMinusMatrix(i,0) - evaluationSurfaceMinusPIA(i,0);
                    controlPointsMinusPIAGPU[i*3+1] = controlPointsMinusPIAGPU[i*3+1] + surfaceMinusMatrix(i,1) - evaluationSurfaceMinusPIA(i,1);
                    controlPointsMinusPIAGPU[i*3+2] = controlPointsMinusPIAGPU[i*3+2] + surfaceMinusMatrix(i,2) - evaluationSurfaceMinusPIA(i,2);
                }
            }
        }

        // compute C k+1
        blocksPerGrid = std::ceil(numberOfPointCloud*1.0/threadsPerBlock.x);
        cudastatus = cudaMemcpy(d_controlPoints, controlPointsPlusPIAGPU, numberControlPointUDirection*numberControlPointVDirection*3*sizeof(float), cudaMemcpyHostToDevice);
        if (cudaSuccess != cudastatus) {
            std::cout << "control points 1 transfer H2D error: " << cudaGetErrorString(cudastatus) << std::endl;
        }
        evaluationGPU<<<blocksPerGrid, threadsPerBlock>>>(d_knot, d_u, d_v, d_controlPoints, d_result, numberOfPointCloud, numberControlPointUDirection, numberControlPointVDirection);

        cudastatus = cudaMemcpy(result, d_result, numberOfPointCloud*3*sizeof(float), cudaMemcpyDeviceToHost);
        if (cudaSuccess != cudastatus) {
            std::cout << "result 1 transfer D2H error: " << cudaGetErrorString(cudastatus) << std::endl;
        }

        for (int i=0;i<numberOfPointCloud;i++) {
            evaluationSurfacePlusPIA(i, 0) = result[i*3];
            evaluationSurfacePlusPIA(i, 1) = result[i*3+1];
            evaluationSurfacePlusPIA(i, 2) = result[i*3+2];
        }

        cudaDeviceSynchronize();
        // if (modelType != "SchwarzP") {
        if (true) {
            cudastatus = cudaMemcpy(d_controlPoints, controlPointsMinusPIAGPU, numberControlPointUDirection*numberControlPointVDirection*3*sizeof(float), cudaMemcpyHostToDevice);
            if (cudaSuccess != cudastatus) {
                std::cout << "control points 2 transfer H2D error: " << cudaGetErrorString(cudastatus) << std::endl;
            }

            evaluationGPU<<<blocksPerGrid, threadsPerBlock>>>(d_knot, d_u, d_v, d_controlPoints, d_result, numberOfPointCloud, numberControlPointUDirection, numberControlPointVDirection);
            
            cudastatus = cudaMemcpy(result, d_result, numberOfPointCloud*3*sizeof(float), cudaMemcpyDeviceToHost);
            if (cudaSuccess != cudastatus) {
                std::cout << "result 2 transfer D2H error: " << cudaGetErrorString(cudastatus) << std::endl;
            }
            if (modelType != "SchwarzP") {
                for (int i=0;i<numberOfPointCloud;i++) {
                    evaluationSurfaceMinusPIA(i, 0) = result[i*3];
                    evaluationSurfaceMinusPIA(i, 1) = result[i*3+1];
                    evaluationSurfaceMinusPIA(i, 2) = result[i*3+2];
                }
            }
        }
        
    }
    for (int i=0;i<numberControlPointUDirection*numberControlPointVDirection;i++) {
        (*controlPointsPlusPIA)(i,0) = controlPointsPlusPIAGPU[i*3];
        (*controlPointsPlusPIA)(i,1) = controlPointsPlusPIAGPU[i*3+1];
        (*controlPointsPlusPIA)(i,2) = controlPointsPlusPIAGPU[i*3+2];
        if (modelType != "SchwarzP") {
            (*controlPointsMinusPIA)(i,0) = controlPointsMinusPIAGPU[i*3];
            (*controlPointsMinusPIA)(i,1) = controlPointsMinusPIAGPU[i*3+1];
            (*controlPointsMinusPIA)(i,2) = controlPointsMinusPIAGPU[i*3+2];
        }
    }

    free(result);
    free(controlPointsPlusPIAGPU);
    free(controlPointsMinusPIAGPU);
    free(uVector);
    free(vVector);

    free(firstOrderEdge1);
    free(firstOrderEdge1_origin);
    free(firstOrderEdge1_second);
    free(firstOrderEdge1_second_origin);
    free(firstOrderEdge2);
    free(firstOrderEdge2_origin);
    free(firstOrderEdge2_second);
    free(firstOrderEdge2_second_origin);
    free(firstOrderEdge3);
    free(firstOrderEdge3_origin);
    free(firstOrderEdge3_second);
    free(firstOrderEdge3_second_origin);
    free(firstOrderEdge4);
    free(firstOrderEdge4_origin);
    free(firstOrderEdge4_second);
    free(firstOrderEdge4_second_origin);
    free(secondOrderEdge1);
    free(secondOrderEdge1_2);
    free(secondOrderEdge1_3);
    free(secondOrderEdge1_origin);
    free(secondOrderEdge1_second_origin);
    free(secondOrderEdge2);
    free(secondOrderEdge2_2);
    free(secondOrderEdge2_3);
    free(secondOrderEdge2_origin);
    free(secondOrderEdge2_second_origin);
    free(secondOrderEdge3);
    free(secondOrderEdge3_2);
    free(secondOrderEdge3_3);
    free(secondOrderEdge3_origin);
    free(secondOrderEdge3_second_origin);
    free(secondOrderEdge4);
    free(secondOrderEdge4_2);
    free(secondOrderEdge4_3);
    free(secondOrderEdge4_origin);
    free(secondOrderEdge4_second_origin);
    cout << "TPMS2STEP > constrained-PIA over" << endl;
}

extern "C"
void memoryAllocation1(float* a, float* b, float* c, float* d, int number) {
    a = (float*)malloc(number*sizeof(float));
    b = (float*)malloc(number*sizeof(float));
    c = (float*)malloc(number*sizeof(float));
    d = (float*)malloc(number*sizeof(float));
    if (!a || !b || !c || !d) {
        cout << "TPMS2STEP > Memory allocation error." << endl;
    }
}

extern "C"
void memoryAllocation2(float* a, float* b, float* c, float* d, float* e, int number) {
    a = (float*)malloc(number*sizeof(float));
    b = (float*)malloc(number*sizeof(float));
    c = (float*)malloc(number*sizeof(float));
    d = (float*)malloc(number*sizeof(float));
    e = (float*)malloc(number*sizeof(float));
}

extern "C"
void myFree1(float* a, float* b, float* c, float* d) {
    if (a) {free(a);}
    if (b) {free(b);}
    if (c) {free(c);}
    if (d) {free(d);}
}

extern "C"
void myFree2(float* a, float* b, float* c, float* d, float* e) {
    if (a) {free(a);}
    if (b) {free(b);}
    if (c) {free(c);}
    if (d) {free(d);}
    if (e) {free(e);}
}

extern "C"
void initialzeCtrlPts(int number, float* ctrlpts, Eigen::MatrixXf& surface) {
    for (int i=0;i<number;i++) {
        ctrlpts[i*3] = surface(i,0);
        ctrlpts[i*3+1] = surface(i,1);
        ctrlpts[i*3+2] = surface(i,2);
    }
}