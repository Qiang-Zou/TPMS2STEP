#include "sample.h"

void sample(std::vector<float>& params) {
    int n = (params.size()+1)/2;
    float length;
    float totalLength = 0;
    std::vector<float> vecLen;
    std::vector<float> vecNormalLen;
    ifstream ifs;
    ifs.open("../data/uv/sampleParams.txt", ios::in);
    for (int i = 0; i < n-1; i++) {
        ifs >> length;
        totalLength += length;
        vecLen.push_back(length);
    }
    float temp = 0;
    for (int i = 0; i < n; i++) {
        vecNormalLen.push_back(temp/totalLength);
        if (i != n-1) {
            temp += vecLen[n-2-i];
        }
    }
    float totalLength2 = 0;
    vecLen.clear();
    std::vector<float> uDirNormalLen;
    for (int i=0;i<n-1;i++) {
        ifs >> length;
        totalLength2 += length;
        vecLen.push_back(length);
    }
    temp = 0;
    for (int i = 0; i < n; i++) {
        uDirNormalLen.push_back(temp/totalLength2);
        if (i != n-1) {
            temp += vecLen[n-2-i];
        }
    }
    ifs.close();
    // for (int i=0;i<params.size();i++) {
    //     cout << "params: " << params[i] << endl;
    // }
    std::vector<float> u1;
    std::vector<float> u2;
    int cursor = 0;
    for (int i=0;i<n-1;i++) {
        while(cursor<=n-1 && vecNormalLen[cursor]<2*params[i]){
            cursor++;
        }
        if (i == 0) {
            u1.push_back(0);
        } else {
            u1.push_back(2*params[cursor-1]+(2*params[i]-vecNormalLen[cursor-1])/(vecNormalLen[cursor]-vecNormalLen[cursor-1])*(2*params[cursor]-2*params[cursor-1]));
        }
    }
    u1.push_back(1);
    cursor = 0;
    for (int i=0;i<n-1;i++) {
        while (cursor<=n-1 && uDirNormalLen[cursor]<2*params[i]) {
            cursor++;
        }
        if (i == 0) {
            u2.push_back(0);
        } else {
            u2.push_back(2*params[cursor-1]+(2*params[i]-uDirNormalLen[cursor-1])/(uDirNormalLen[cursor]-uDirNormalLen[cursor-1])*(2*params[cursor]-2*params[cursor-1]));
        }
    }
    u2.push_back(1);

    for (int i=0;i<n;i++) {
        u1[i] = u1[i]*0.5;
    }
    for (int i=0;i<n-1;i++) {
        u1.push_back(1.0-u1[n-2-i]);
    }
    ofstream ofs;
    ofs.open("../data/uv/uv_Gyroid.off", ios::out);
    ofs << "OFF\n" << std::to_string((int)(n*n)) << " 0 0" << endl;
    float* endVerticesX = (float*)malloc((2*n-1) * sizeof(float));
    float* endVerticesY = (float*)malloc((2*n-1) * sizeof(float));
    for (int i = 0; i < n; i++) {
        if (i == 1) {
            float sinTheta = (sqrt(6)-sqrt(2))/4*(1-u1[i]*0.5/0.5);
            float cosTheta = sqrt(1 - sinTheta * sinTheta);
            endVerticesX[i] = cosTheta + sinTheta - 1 / sqrt(2);
            endVerticesY[i] = cosTheta - sinTheta - 1 / sqrt(2);
        } else if(i == 2) {
            float sinTheta = (sqrt(6)-sqrt(2))/4*(1-u1[i]*0.8/0.5);
            float cosTheta = sqrt(1 - sinTheta * sinTheta);
            endVerticesX[i] = cosTheta + sinTheta - 1 / sqrt(2);
            endVerticesY[i] = cosTheta - sinTheta - 1 / sqrt(2);
        } else {
            float sinTheta = (sqrt(6)-sqrt(2))/4*(1-u1[i]/0.5);
            float cosTheta = sqrt(1 - sinTheta * sinTheta);
            endVerticesX[i] = cosTheta + sinTheta - 1 / sqrt(2);
            endVerticesY[i] = cosTheta - sinTheta - 1 / sqrt(2);
        }
    }
    for (int i = n; i < 2*n-1; i++) {
        if (i == 2*n-3) {
            float sinTheta = (sqrt(6)-sqrt(2))/4*(1-(1-u1[i])*0.5/0.5);
            float cosTheta = sqrt(1 - sinTheta * sinTheta);
            endVerticesX[i] = cosTheta - sinTheta - 1 / sqrt(2);
            endVerticesY[i] = cosTheta + sinTheta - 1 / sqrt(2);
        } else if (i == 2*n-4) {
            float sinTheta = (sqrt(6)-sqrt(2))/4*(1-(1-u1[i])*0.8/0.5);
            float cosTheta = sqrt(1 - sinTheta * sinTheta);
            endVerticesX[i] = cosTheta - sinTheta - 1 / sqrt(2);
            endVerticesY[i] = cosTheta + sinTheta - 1 / sqrt(2);
        } else {
            float sinTheta = (sqrt(6)-sqrt(2))/4*(1-(1-u1[i])/0.5);
            float cosTheta = sqrt(1 - sinTheta * sinTheta);
            endVerticesX[i] = cosTheta - sinTheta - 1 / sqrt(2);
            endVerticesY[i] = cosTheta + sinTheta - 1 / sqrt(2);
        }
    }
    for (int i=0;i<2*n-1;i++) {
        ofs << endVerticesX[i] << " " << endVerticesY[i] << " 0" << endl;
        // cout << endVerticesX[i] << " " << endVerticesY[i] << endl;
    }
    float *midVertices = (float*)malloc(n*sizeof(float));
    for (int i=0;i<n-1;i++) {
        midVertices[i] = (sqrt(2)-1)*(1-u2[i])/sqrt(2);
        // cout << "midVertices: " << midVertices[i] << endl;
    }
    float le = (sqrt(3)-1)/sqrt(2);
    for (int i=1;i<n;i++) {
        float paramU = 1-2*u1[i];
        std::vector<float> tmpx;
        std::vector<float> tmpy;
        for (int j=i;j<n-1;j++) {
            float paramV = 1-2*u1[j];
            float alpha = midVertices[i]/(midVertices[i]-le*paramU);
            float x0 = (alpha*le*paramU-le*paramV)/(alpha-1);
            float y0 = x0-le*paramV;
            ofs << x0 << " " << y0 << " 0" << endl;
            tmpx.push_back(x0);
            tmpy.push_back(y0);
        }
        if (i!=n-1) {
            ofs << midVertices[i] << " " << midVertices[i] << " 0" << endl;
            for (int k=n-i-2;k>=0;k--) {
                ofs << tmpy[k] << " " << tmpx[k] << " 0" << endl;
            }
        }
    }
    ofs << "0 0 0" << endl;
    ofs.close();
    free(endVerticesX);
    free(endVerticesY);
}

void sample2(std::vector<float>& params) {
    int n = (params.size()+1)/2;
    float delta = sqrt(2)-1;
    float theta = M_PI/6/(2*n-2);
    ofstream ofs;
    ofs.open("../data/uv/uv_Gyroid.off", ios::out);
    ofs << "OFF\n" << std::to_string((int)(n*n)) << " 0 0" << endl;
    for (int i=0;i<n;i++) {
        float ri = 1.0/(sqrt(2)*sin(theta*params[i]+M_PI/6));
        for (int j=i;j<2*n-1-i;j++) {
            float thetaj = M_PI/6+theta*j;
            float x0 = ri*cos(thetaj)-1/sqrt(2);
            float y0 = ri*sin(thetaj)-1/sqrt(2);
            ofs << x0 << " " << y0 << " 0" << endl;
        }
    }
    ofs.close();
}

int findSpan2(int degree, int n, float t, float* knot) {
    if (t == knot[n+1]) {
        return n;
    }
    int low = degree;
    int high = n + 1;
    int mid;

    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (t < knot[mid]-1e-5)
            high = mid;
        else
            low = mid;
    }
    return low;
}

void sampleDiamond(std::vector<float>& params) {
    int n = (params.size()+1)/2;
    float delta = sqrt(2)-1;
    float theta = M_PI/12;
    ofstream ofs;
    ofs.open("../data/uv/uv_Diamond.off", ios::out);
    ofs << "OFF\n" << std::to_string((int)(n*(n+1)/2)) << " 0 0" << endl;
    for (int i=0;i<n;i++) {
        float ri = sqrt(((sqrt(3)-1)/sqrt(2)*(1-2*params[i])+1.0/sqrt(2))*((sqrt(3)-1)/sqrt(2)*(1-2*params[i])+1.0/sqrt(2))+0.5);
        for (int j=0;j<n-i;j++) {
            float  rj = sqrt(((sqrt(3)-1)/sqrt(2)*(1-2*params[j+i])+1.0/sqrt(2))*((sqrt(3)-1)/sqrt(2)*(1-2*params[j+i])+1.0/sqrt(2))+0.5);
            float thetaj = asin(1.0/(sqrt(2)*rj));
            // cout << "ri: " << ri << " thetaj: " << thetaj << " " << asin(1.0/sqrt(2)*ri) << " " << theta*2*params[j] << endl;
            float x0 = ri*cos(thetaj)-1/sqrt(2);
            float y0 = ri*sin(thetaj)-1/sqrt(2);
            // cout << x0 << " " << y0 << " 0" << endl;
            ofs << x0 << " " << y0 << " 0" << endl;
        }
    }
    ofs.close();
}

void sampleSchwarzP(std::vector<float>& params) {
    int n = (params.size()+1)/2;
    float delta = sqrt(2)-1;
    float theta = M_PI/12;
    ofstream ofs;
    ofs.open("../data/uv/uv_SchwarzP.off", ios::out);
    ofs << "OFF\n" << std::to_string((int)(n*(n+1)/2)) << " 0 0" << endl;
    for (int i=0;i<n;i++) {
        float ri = sqrt(((sqrt(3)-1)/sqrt(2)*(1-2*params[i])+1.0/sqrt(2))*((sqrt(3)-1)/sqrt(2)*(1-2*params[i])+1.0/sqrt(2))+0.5);
        for (int j=0;j<n-i;j++) {
            float  rj = sqrt(((sqrt(3)-1)/sqrt(2)*(1-2*params[j+i])+1.0/sqrt(2))*((sqrt(3)-1)/sqrt(2)*(1-2*params[j+i])+1.0/sqrt(2))+0.5);
            float thetaj = asin(1.0/(sqrt(2)*rj));
            float x0 = ri*cos(thetaj)-1/sqrt(2);
            float y0 = ri*sin(thetaj)-1/sqrt(2);
            ofs << x0 << " " << y0 << " 0" << endl;
            
        }
    }
    ofs.close();
}