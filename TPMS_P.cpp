#include "TPMS_P.h"

TPMS_P::TPMS_P()
{
}

TPMS_P::~TPMS_P()
{
}

bool TPMS_P::constructSolidModel() {
    this->cellScaleX = 1.0;
    this->cellScaleY = 1.0;
    this->cellScaleZ = 1.0;
    this->cellSizeX = 2.0 * this->cellScaleX;
    this->cellSizeY = 2.0 * this->cellScaleY;
    this->cellSizeZ = 2.0 * this->cellScaleZ;
    this->maxIteration = 1;
    this->setOffsetValue(0.1, 0);
    this->setCellNumbers(1, 1, 1);
    this->setTolerance(0.01);
    float epsilon_1 = 0.5*this->getTolerance();
    float epsilon_2 = 0.5*this->getTolerance();

    // define some params
    ifstream ifs;
    ofstream ofs;
    int UDegree = 3;
    int VDegree = 3;
    int NumRows;
    int NumCols;
    Eigen::VectorXf U;
    Eigen::VectorXf V;
    int numberOfPointCloud = 0;
    Eigen::MatrixXf surfacePlusMatrix;
    Eigen::MatrixXf surfaceMinusMatrix;
    Eigen::MatrixXf controlPointsPlusPIA;
    Eigen::MatrixXf controlPointsMinusPIA;
    this->minControlPointNumberPerEdge = 10;
    float* knot = NULL;

    // estimate the sample density 
    float offset = std::max(std::fabs(this->getOffsetValue1()), std::fabs(this->getOffsetValue2()));
    float m1 = 2.0 + 4.0 * std::fabs(offset);
    float m2 = 2.0 + 4.0 * std::fabs(offset);
    float m3 = 2.0 + 4.0 * std::fabs(offset);
    float L = sqrt((8 * epsilon_1) / (m1 + 2 * m2 + m3));
    int n = std::floor((std::sqrt(2) - 1) / L) + 1;
    if (n<this->minControlPointNumberPerEdge) {
        n=this->minControlPointNumberPerEdge;
    }
    std::vector<float> paramRuler;
    vector<float> knotTemp;
    errorControlledSampling(knotTemp, knot, &NumRows, &NumCols, &numberOfPointCloud, surfacePlusMatrix, surfaceMinusMatrix, U, V, controlPointsPlusPIA, controlPointsMinusPIA, n, epsilon_2, paramRuler);
    cout << "TPMS2STEP > Error Control ok." << endl;
    if (knot) {
        free(knot);
    }
    knot = (float*)malloc(knotTemp.size()*sizeof(float));
    for (int i=0;i<knotTemp.size();i++) {
        knot[i] = knotTemp[i];
    }

    // define knot vectors
    TColStd_Array1OfReal* Uknot = new TColStd_Array1OfReal(1, (paramRuler.size()+1)/2-2);
    TColStd_Array1OfReal* Vknot = new TColStd_Array1OfReal(1, (paramRuler.size()+1)/2-2);

    for (int i=1;i<=(paramRuler.size()+1)/2-2;i++) {
        (*Uknot).SetValue(i,(i-1)*1.0/((paramRuler.size()+1)/2-3));
        (*Vknot).SetValue(i,(i-1)*1.0/((paramRuler.size()+1)/2-3));
    }
    int tempLength = (paramRuler.size()+1)/2;
    TColStd_Array1OfInteger* multiplicitiesU = new TColStd_Array1OfInteger(1, tempLength-2);
    TColStd_Array1OfInteger* multiplicitiesV = new TColStd_Array1OfInteger(1, tempLength-2);
    for (int i=2;i<=tempLength-3;i++) {
        (*multiplicitiesU).SetValue(i,1);
        (*multiplicitiesV).SetValue(i,1);
    }
    (*multiplicitiesU).SetValue(1, 4);
    (*multiplicitiesU).SetValue(tempLength-2, 4);
    (*multiplicitiesV).SetValue(1, 4);
    (*multiplicitiesV).SetValue(tempLength-2, 4);

    TColgp_Array2OfPnt* ctrlpts_1 = new TColgp_Array2OfPnt(1, NumRows, 1, NumCols);
    TColgp_Array2OfPnt* ctrlpts_2 = new TColgp_Array2OfPnt(1, NumRows, 1, NumCols);

    // assembling

#ifdef USE_QUAD_HEX_MESH
    int numberOfQuadMesh = 8;
    std::vector<float> quadMeshPoints_1;
    std::vector<float> quadMeshPoints_2;
    std::vector<float> quadMeshPointsTotal;
#endif

    const double ksi = 1.0 / sqrt(2.0);
    std::vector<math_Matrix> F;
    math_Matrix F1(1, 3, 1, 3, 0);
    math_Matrix F2(1, 3, 1, 3, 0);
    math_Matrix F3(1, 3, 1, 3, 0);
    math_Matrix F4(1, 3, 1, 3, 0);
    math_Matrix F5(1, 3, 1, 3, 0);
    math_Matrix F6(1, 3, 1, 3, 0);
    F1(1, 1) = (-1) * ksi; F1(1, 2) = (-1) * ksi; F1(1, 3) = 0;
    F1(2, 1) = 0; F1(2, 2) = 0; F1(2, 3) = 1;
    F1(3, 1) = (-1) * ksi; F1(3, 2) = ksi; F1(3, 3) = 0;

    F2(1, 1) = (-1) * ksi; F2(1, 2) = (-1) * ksi; F2(1, 3) = 0;
    F2(2, 1) = (-1) * ksi; F2(2, 2) = ksi; F2(2, 3) = 0;
    F2(3, 1) = 0; F2(3, 2) = 0; F2(3, 3) = 1;

    F3(1, 1) = 0; F3(1, 2) = 0; F3(1, 3) = 1;
    F3(2, 1) = (-1) * ksi; F3(2, 2) = (-1) * ksi; F3(2, 3) = 0;
    F3(3, 1) = (-1) * ksi; F3(3, 2) = ksi; F3(3, 3) = 0;

    F4(1, 1) = (-1) * ksi; F4(1, 2) = ksi; F4(1, 3) = 0;
    F4(2, 1) = 0; F4(2, 2) = 0; F4(2, 3) = 1;
    F4(3, 1) = (-1) * ksi; F4(3, 2) = (-1) * ksi; F4(3, 3) = 0;

    F5(1, 1) = (-1) * ksi; F5(1, 2) = ksi; F5(1, 3) = 0;
    F5(2, 1) = (-1) * ksi; F5(2, 2) = (-1) * ksi; F5(2, 3) = 0;
    F5(3, 1) = 0; F5(3, 2) = 0; F5(3, 3) = 1;

    F6(1, 1) = 0; F6(1, 2) = 0; F6(1, 3) = 1;
    F6(2, 1) = (-1) * ksi; F6(2, 2) = ksi; F6(2, 3) = 0;
    F6(3, 1) = (-1) * ksi; F6(3, 2) = (-1) * ksi; F6(3, 3) = 0;

    F.push_back(F1);
    F.push_back(F2);
    F.push_back(F3);
    F.push_back(F4);
    F.push_back(F5);
    F.push_back(F6);

    std::vector<math_Vector> eta_f;
    math_Vector eta_f1(1, 3, 0);
    math_Vector eta_f2(1, 3, 0);
    math_Vector eta_f3(1, 3, 0);
    math_Vector eta_f4(1, 3, 0);
    math_Vector eta_f5(1, 3, 0);
    math_Vector eta_f6(1, 3, 0);

    eta_f1(1) = 0; eta_f1(2) = 0.5; eta_f1(3) = 1;
    eta_f2(1) = 0; eta_f2(2) = 1; eta_f2(3) = 0.5;
    eta_f3(1) = 0.5; eta_f3(2) = 0; eta_f3(3) = 1;
    eta_f4(1) = 1; eta_f4(2) = 0.5; eta_f4(3) = 0;
    eta_f5(1) = 1; eta_f5(2) = 0; eta_f5(3) = 0.5;
    eta_f6(1) = 0.5; eta_f6(2) = 1; eta_f6(3) = 0;

    eta_f.push_back(eta_f1);
    eta_f.push_back(eta_f2);
    eta_f.push_back(eta_f3);
    eta_f.push_back(eta_f4);
    eta_f.push_back(eta_f5);
    eta_f.push_back(eta_f6);

    std::vector<math_Matrix> Q;
    std::vector<math_Vector> eta;
    math_Matrix Q1(1, 3, 1, 3, 0);
    math_Matrix Q2(1, 3, 1, 3, 0);
    math_Matrix Q3(1, 3, 1, 3, 0);
    math_Matrix Q4(1, 3, 1, 3, 0);
    math_Matrix Q5(1, 3, 1, 3, 0);
    math_Matrix Q6(1, 3, 1, 3, 0);
    math_Matrix Q7(1, 3, 1, 3, 0);
    math_Matrix Q8(1, 3, 1, 3, 0);
    math_Vector eta1(1, 3, 0);
    math_Vector eta2(1, 3, 0);
    math_Vector eta3(1, 3, 0);
    math_Vector eta4(1, 3, 0);
    math_Vector eta5(1, 3, 0);
    math_Vector eta6(1, 3, 0);
    math_Vector eta7(1, 3, 0);
    math_Vector eta8(1, 3, 0);

    Q1(1, 1) = 1; Q1(2, 2) = 1; Q1(3, 3) = 1;
    Q2(1, 1) = -1; Q2(2, 2) = 1; Q2(3, 3) = 1;
    Q3(1, 1) = 1; Q3(2, 2) = -1; Q3(3, 3) = 1;
    Q4(1, 1) = -1; Q4(2, 2) = -1; Q4(3, 3) = 1;
    Q5(1, 1) = 1; Q5(2, 2) = 1; Q5(3, 3) = -1;
    Q6(1, 1) = -1; Q6(2, 2) = 1; Q6(3, 3) = -1;
    Q7(1, 1) = 1; Q7(2, 2) = -1; Q7(3, 3) = -1;
    Q8(1, 1) = -1; Q8(2, 2) = -1; Q8(3, 3) = -1;

    Q.push_back(Q1);
    Q.push_back(Q2);
    Q.push_back(Q3);
    Q.push_back(Q4);
    Q.push_back(Q5);
    Q.push_back(Q6);
    Q.push_back(Q7);
    Q.push_back(Q8);

    eta.push_back(eta1);
    eta.push_back(eta2);
    eta.push_back(eta3);
    eta.push_back(eta4);
    eta.push_back(eta5);
    eta.push_back(eta6);
    eta.push_back(eta7);
    eta.push_back(eta8);

    math_Matrix surface_0(1, (*ctrlpts_1).Size(), 1, 3);
    std::vector<TColgp_Array2OfPnt> controlPointsVector;
    std::vector<math_Matrix> surfaceMatrixVector;
    GeomAPI_PointsToBSplineSurface bsplineSurface;
    std::vector<TColgp_Array2OfPnt> controlPointsTotalVector;
    std::vector<TColgp_Array1OfPnt> curveControlPointsVector;

    for (int i=1;i<=NumRows;i++) {
        for (int j=1;j<=NumCols;j++) {
            (*ctrlpts_1).SetValue(i,j,gp_Pnt(controlPointsPlusPIA((i-1)*NumCols+(j-1), 0), controlPointsPlusPIA((i-1)*NumCols+(j-1), 1), controlPointsPlusPIA((i-1)*NumCols+(j-1), 2)));
            (*ctrlpts_2).SetValue(i,j,gp_Pnt(controlPointsMinusPIA((i-1)*NumCols+(j-1), 0), controlPointsMinusPIA((i-1)*NumCols+(j-1), 1), controlPointsMinusPIA((i-1)*NumCols+(j-1), 2)));
            surface_0((i-1)*NumCols+j, 1) = controlPointsPlusPIA((i-1)*NumCols+(j-1), 0);
            surface_0((i-1)*NumCols+j, 2) = controlPointsPlusPIA((i-1)*NumCols+(j-1), 1);
            surface_0((i-1)*NumCols+j, 3) = controlPointsPlusPIA((i-1)*NumCols+(j-1), 2);
        }
    }
    int controlPointsColLength = (*ctrlpts_1).ColLength();
    int controlPointsRowLength = (*ctrlpts_1).RowLength();

#ifdef USE_QUAD_HEX_MESH
    Handle(Geom_BSplineSurface) firstSplineSurface = new Geom_BSplineSurface(
        *ctrlpts_1,
        *Uknot,
        *Vknot,
        *multiplicitiesU,
        *multiplicitiesV,
        3,
        3
    );
    Handle(Geom_BSplineSurface) secondSplineSurface = new Geom_BSplineSurface(
        *ctrlpts_2,
        *Uknot,
        *Vknot,
        *multiplicitiesU,
        *multiplicitiesV,
        3,
        3
    );
    extractMesh(firstSplineSurface, quadMeshPoints_1, numberOfQuadMesh);
    extractMesh(secondSplineSurface, quadMeshPoints_2, numberOfQuadMesh);

    math_Matrix quadMeshMatrix_1(1, numberOfQuadMesh*numberOfQuadMesh, 1, 3);
    math_Matrix quadMeshMatrix_2(1, numberOfQuadMesh*numberOfQuadMesh, 1, 3);
    std::vector<math_Matrix> quadMeshVector_1;

    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,1) = quadMeshPoints_1[(i*numberOfQuadMesh+j)*3];
            quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,2) = quadMeshPoints_1[(i*numberOfQuadMesh+j)*3+1];
            quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,3) = quadMeshPoints_1[(i*numberOfQuadMesh+j)*3+2];
            quadMeshMatrix_2(i*numberOfQuadMesh+j+1 ,1) = quadMeshPoints_2[(i*numberOfQuadMesh+j)*3];
            quadMeshMatrix_2(i*numberOfQuadMesh+j+1 ,2) = quadMeshPoints_2[(i*numberOfQuadMesh+j)*3+1];
            quadMeshMatrix_2(i*numberOfQuadMesh+j+1 ,3) = quadMeshPoints_2[(i*numberOfQuadMesh+j)*3+2];
        }
    }
#endif

    for (int k = 0; k < 6; k++) {
        math_Matrix sampleSurface(1, (*ctrlpts_1).Size(), 1, 3);
        math_Matrix tempS = F[k] * surface_0.Transposed();
        math_Matrix tempSS = tempS.Transposed();
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                sampleSurface((i - 1) * NumCols + j, 1) = (tempSS((i - 1) * NumCols + j, 1) + eta_f[k](1));
                sampleSurface((i - 1) * NumCols + j, 2) = (tempSS((i - 1) * NumCols + j, 2) + eta_f[k](2));
                sampleSurface((i - 1) * NumCols + j, 3) = (tempSS((i - 1) * NumCols + j, 3) + eta_f[k](3));
            }
        }

        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                (*ctrlpts_1).SetValue(i, j, gp_Pnt(sampleSurface((i - 1) * controlPointsRowLength + j, 1), sampleSurface((i - 1) * controlPointsRowLength + j, 2), sampleSurface((i - 1) * controlPointsRowLength + j, 3)));
            }
        }
        controlPointsVector.push_back((*ctrlpts_1));
        surfaceMatrixVector.push_back(sampleSurface);
        TColgp_Array1OfPnt curveCpts(1, controlPointsColLength);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            curveCpts.SetValue(i, (*ctrlpts_1).Value(controlPointsColLength, i));
        }
        curveControlPointsVector.push_back(curveCpts);

#ifdef USE_QUAD_HEX_MESH
        math_Matrix quadTemp = F[k] * quadMeshMatrix_1.Transposed();
        math_Matrix quadTemp2 = quadTemp.Transposed();
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadTemp2(i*numberOfQuadMesh+j+1, 1) += eta_f[k](1);
                quadTemp2(i*numberOfQuadMesh+j+1, 2) += eta_f[k](2);
                quadTemp2(i*numberOfQuadMesh+j+1, 3) += eta_f[k](3);
            }
        }
        quadMeshVector_1.push_back(quadTemp2);
#endif

    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 8; j++) {
            math_Matrix surfaceTemp(1, controlPointsColLength * controlPointsRowLength, 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3);
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(controlPointsTemp);

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp = Q[j] * quadMeshVector_1[i].Transposed();
            math_Matrix quadTemp2 = quadTemp.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    controlPointsVector.clear();
    surfaceMatrixVector.clear();
    
#ifdef USE_QUAD_HEX_MESH
    quadMeshVector_1.clear();
    quadMeshPoints_1.clear();
#endif

    // read the surface data
    readSemiFundamentalPatch(surfacePlusMatrix, surfaceMinusMatrix, &numberOfPointCloud, this->getOffsetValue2());

    constrainedPIAGPU(knot, 3, numberOfPointCloud, NumRows, NumCols, surfacePlusMatrix, surfaceMinusMatrix, U, V,
    &controlPointsPlusPIA, &controlPointsMinusPIA, this->getOffsetValue2(), this->modelType);
    for (int i=1;i<=NumRows;i++) {
        for (int j=1;j<=NumCols;j++) {
            (*ctrlpts_1).SetValue(i,j,gp_Pnt(controlPointsPlusPIA((i-1)*NumCols+(j-1), 0), controlPointsPlusPIA((i-1)*NumCols+(j-1), 1), controlPointsPlusPIA((i-1)*NumCols+(j-1), 2)));
            surface_0((i-1)*NumCols+j, 1) = controlPointsPlusPIA((i-1)*NumCols+(j-1), 0);
            surface_0((i-1)*NumCols+j, 2) = controlPointsPlusPIA((i-1)*NumCols+(j-1), 1);
            surface_0((i-1)*NumCols+j, 3) = controlPointsPlusPIA((i-1)*NumCols+(j-1), 2);
        }
    }
    controlPointsColLength = (*ctrlpts_1).ColLength();
    controlPointsRowLength = (*ctrlpts_1).RowLength();

#ifdef USE_QUAD_HEX_MESH
    Handle(Geom_BSplineSurface) secondBSplineSurface = new Geom_BSplineSurface(
        *ctrlpts_1,
        *Uknot,
        *Vknot,
        *multiplicitiesU,
        *multiplicitiesV,
        3,
        3
    );
    extractMesh(secondBSplineSurface, quadMeshPoints_1, numberOfQuadMesh);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,1) = quadMeshPoints_1[(i*numberOfQuadMesh+j)*3];
            quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,2) = quadMeshPoints_1[(i*numberOfQuadMesh+j)*3+1];
            quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,3) = quadMeshPoints_1[(i*numberOfQuadMesh+j)*3+2];
        }
    }
#endif

    for (int k = 0; k < 6; k++) {
        math_Matrix sampleSurface(1, (*ctrlpts_1).Size(), 1, 3);
        math_Matrix tempS = F[k] * surface_0.Transposed();
        math_Matrix tempSS = tempS.Transposed();
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                sampleSurface((i - 1) * controlPointsColLength + j, 1) = (tempSS((i - 1) * controlPointsColLength + j, 1) + eta_f[k](1));
                sampleSurface((i - 1) * controlPointsColLength + j, 2) = (tempSS((i - 1) * controlPointsColLength + j, 2) + eta_f[k](2));
                sampleSurface((i - 1) * controlPointsColLength + j, 3) = (tempSS((i - 1) * controlPointsColLength + j, 3) + eta_f[k](3));
            }
        }
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                (*ctrlpts_1).SetValue(i, j, gp_Pnt(sampleSurface((i - 1) * controlPointsRowLength + j, 1), sampleSurface((i - 1) * controlPointsRowLength + j, 2), sampleSurface((i - 1) * controlPointsRowLength + j, 3)));
            }
        }
        controlPointsVector.push_back((*ctrlpts_1));
        surfaceMatrixVector.push_back(sampleSurface);
        TColgp_Array1OfPnt curveCpts(1, controlPointsColLength);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            curveCpts.SetValue(i, (*ctrlpts_1).Value(controlPointsColLength, i));
        }
        curveControlPointsVector.push_back(curveCpts);

#ifdef USE_QUAD_HEX_MESH
        math_Matrix quadTemp = F[k] * quadMeshMatrix_1.Transposed();
        math_Matrix quadTemp2 = quadTemp.Transposed();
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadTemp2(i*numberOfQuadMesh+j+1, 1) += eta_f[k](1);
                quadTemp2(i*numberOfQuadMesh+j+1, 2) += eta_f[k](2);
                quadTemp2(i*numberOfQuadMesh+j+1, 3) += eta_f[k](3);
            }
        }
        quadMeshVector_1.push_back(quadTemp2);
#endif

    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 8; j++) {
            math_Matrix surfaceTemp(1, controlPointsColLength * controlPointsRowLength, 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3);
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(controlPointsTemp);

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp = Q[j] * quadMeshVector_1[i].Transposed();
            math_Matrix quadTemp2 = quadTemp.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    std::vector<Handle(Geom_BSplineSurface)> cellSurfaceVector;
    std::vector<Handle(Geom_BSplineCurve)> bSplineCurveVector;

    for (int i = 0; i < 96; i++) {
        for (int j=1;j<=controlPointsTotalVector[i].ColLength();j++) {
            for (int k=1;k<=controlPointsTotalVector[i].RowLength();k++) {
                controlPointsTotalVector[i].SetValue(j,k,gp_Pnt(controlPointsTotalVector[i].Value(j,k).X()*cellScaleX, controlPointsTotalVector[i].Value(j,k).Y()*cellScaleY, controlPointsTotalVector[i].Value(j,k).Z()*cellScaleZ));
            }
        }
        Handle(Geom_BSplineSurface) bsplineTemp = new Geom_BSplineSurface(
            controlPointsTotalVector[i],
            (*Uknot),
            (*Vknot),
            (*multiplicitiesU),
            (*multiplicitiesV),
            3,
            3
        );
        cellSurfaceVector.push_back(bsplineTemp);
    }
    for (int i = 0; i < 12; i++) {
        for (int j=1;j<=curveControlPointsVector[i].Size();j++) {
            curveControlPointsVector[i].SetValue(j, gp_Pnt(curveControlPointsVector[i].Value(j).X()*cellScaleX, curveControlPointsVector[i].Value(j).Y()*cellScaleY, curveControlPointsVector[i].Value(j).Z()*cellScaleZ));
        }
        Handle(Geom_BSplineCurve) bsplineTemp = new Geom_BSplineCurve(
            curveControlPointsVector[i],
            (*Uknot),
            (*multiplicitiesU),
            3
        );
        bSplineCurveVector.push_back(bsplineTemp);
    }

    std::vector<BRepBuilderAPI_MakeFace> faceMakerVector;
    std::vector<TopoDS_Face> topoFaceVector;
    int cellFaceNumber = 96;
    int internalFaceNumber = cellFaceNumber * this->getCellNumberX() * this->getCellNumberY() * this->getCellNumberZ();
    int coverNumber = 8 * 2 * (this->getCellNumberX() * this->getCellNumberY() + this->getCellNumberX() * this->getCellNumberZ() + this->getCellNumberY() * this->getCellNumberZ());
    int faceNumber = internalFaceNumber + coverNumber;
    std::vector<Handle(Geom_BSplineSurface)> faceVector;
    std::vector<Handle(Geom_BSplineCurve)> newCurve;
    std::vector<Geom_BSplineSurface> basicSurfaceElementVector;

    buildCurvesSchwarzP(newCurve, bSplineCurveVector, basicSurfaceElementVector);

    // add the gap faces
    for (int i = 0; i < 8; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[i].Copy());
        for (int i = 0; i < this->getCellNumberX(); i++) {        
            for (int j = 0; j < this->getCellNumberY(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(i * cellSizeX, j * cellSizeY, 0)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, 0, this->getCellNumberZ() * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }
    for (int i = 8; i < 16; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[i].Copy());
        for (int i = 0; i < this->getCellNumberY(); i++) {
            for (int j = 0; j < this->getCellNumberZ(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(0, i * cellSizeY, j * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(this->getCellNumberX() * cellSizeX, 0, 0)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }
    for (int i = 16; i < 24; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[i].Copy());
        for (int i = 0; i < this->getCellNumberZ(); i++) {
            for (int j = 0; j < this->getCellNumberX(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(j * cellSizeX, 0, i * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, this->getCellNumberY() * cellSizeY, 0)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }

#ifdef USE_QUAD_HEX_MESH
    std::vector<float> quadMeshAll;
    for (int i = 0; i < this->getCellNumberX(); i++) {
        for (int j = 0; j < this->getCellNumberY(); j++) {
            for (int k = 0; k < this->getCellNumberZ(); k++) {
                for (int m = 0; m < numberOfQuadMesh*numberOfQuadMesh*96; m++) {
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3]+i*cellSizeX);
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3+1]+j*cellSizeY);
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3+2]+k*cellSizeZ);
                }
            }
        }
    }
    writeHexahedralMesh(quadMeshAll, 96*this->getCellNumberX()*this->getCellNumberY()*this->getCellNumberZ(), 48, numberOfQuadMesh, "schwarzp_hex.off");

    std::vector<float> quadMeshEdge;
    for (int i=0;i<faceVector.size();i++) {
        quadMeshEdge.clear();
        extractMesh(faceVector[i], quadMeshEdge, numberOfQuadMesh);
        for (int j = 0; j < numberOfQuadMesh*numberOfQuadMesh*3; j++) {
            quadMeshAll.push_back(quadMeshEdge[j]);
        }
    }

    writeMesh(quadMeshAll, 96*this->getCellNumberX()*this->getCellNumberY()*this->getCellNumberZ()+faceVector.size(), numberOfQuadMesh, "schwarzp.off");
#endif

    for (int i = 0; i < this->getCellNumberX(); i++) {
        for (int j = 0; j < this->getCellNumberY(); j++) {
            for (int k = 0; k < this->getCellNumberZ(); k++) {
                for (int m = 0; m < cellSurfaceVector.size(); m++) {
                    Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(cellSurfaceVector[m]->Translated(gp_Vec(i * cellSizeX, j * cellSizeY, k * cellSizeZ)));
                    faceVector.push_back(tempTranslatedFace);
                }
            }
        }
    }

    float scaleFactor = 0;
    for (int i = 0; i < faceNumber; i++) {
        /************    scale   **************/
        for (int j=1;j<=faceVector[i]->NbUPoles();j++) {
            for (int k=1;k<=faceVector[i]->NbVPoles();k++) {
                faceVector[i]->SetPole(j,k,gp_Pnt((faceVector[i]->Pole(j,k)).X(),
                (faceVector[i]->Pole(j,k)).Y()*(1+scaleFactor), (faceVector[i]->Pole(j,k)).Z()));
            }
        }
        /**************************/
        BRepBuilderAPI_MakeFace faceMakerTemp(faceVector[i], 1.0e-6);
        faceMakerVector.push_back(faceMakerTemp);
    }

    bool flag = 1;
    for (int i = 0; i < faceNumber; i++) {
        if (!faceMakerVector[i].IsDone()) {
            flag = 0;
            break;
        }
    }
    if (flag) {
        for (int i = 0; i < faceNumber; i++) {
            topoFaceVector.push_back(faceMakerVector[i].Face());
        }
        TopoDS_Builder builder;
        TopoDS_Shell tpmsShell;
        tpmsSolidModel = new TopoDS_Solid();
        builder.MakeShell(tpmsShell);
        for (int i = 0; i < faceNumber; i++) {
            builder.Add(tpmsShell, topoFaceVector[i]);
        }
        BRepBuilderAPI_Sewing aSewingTool;
        aSewingTool.Init(1e-2, Standard_True);
        aSewingTool.Load(tpmsShell);
        aSewingTool.Perform();
        builder.MakeSolid(*tpmsSolidModel);
        builder.Add(*tpmsSolidModel, TopoDS::Shell(aSewingTool.SewedShape()));
        BRepCheck_Analyzer analyzer(*tpmsSolidModel);
        if (analyzer.IsValid()) {
            cout << "TPMS2STEP > topo solid valid" << endl;
        }
        else {
            cout << "TPMS2STEP > topo solid error" << endl;
        }

        Handle(TopTools_HSequenceOfShape) aHSequenceOfShape = new TopTools_HSequenceOfShape;
        aHSequenceOfShape->Clear();
        aHSequenceOfShape->Append(*tpmsSolidModel);
        Standard_CString fileName = outputSTEPFileName.c_str();
        // save file
        IFSelect_ReturnStatus status;
        STEPControl_StepModelType aValue = STEPControl_ManifoldSolidBrep;
        STEPControl_Writer aWriter;

        for (int i = 1; i <= aHSequenceOfShape->Length(); i++) {
            status = aWriter.Transfer(aHSequenceOfShape->Value(i), aValue);
            if (status != IFSelect_RetDone) {
                std::cout << "TPMS2STEP > Write to STEP files error." << endl;
                return 1;
            }
        }
        // write data to disk
        status = aWriter.Write(fileName);
    } else {
        std::cerr << "TPMS2STEP > Failed to create SchwarzP TPMS solid model." << std::endl;
        return false;
    }
    cout << "TPMS2STEP > SchwarzP TPMS solid model generated successfully." << endl;
    return true;
}

void TPMS_P::errorControlledSampling(vector<float>& knotTemp, float* knot, int* NumRows, int* NumCols, int* numberOfPointCloud, Eigen::MatrixXf& surfacePlusMatrix,
Eigen::MatrixXf& surfaceMinusMatrix, Eigen::VectorXf& U, Eigen::VectorXf& V, Eigen::MatrixXf& controlPointsPlusPIA, Eigen::MatrixXf& controlPointsMinusPIA,
int n, float epsilon_2, std::vector<float>& paramRuler) {
    bool flagErrorControl = true;
    int iterationCount = 0;
    string tempString;
    float a, b, c;
    ofstream ofs;
    int gridLength;
    float lGrid;
    float maxNURBSerror;
    std::vector<int> refineTimes;
    for (int i=0;i<2*n;i++) {
        refineTimes.push_back(0);
    }
    // sample grid
    paramRuler.clear();
    for (int i=0;i<2*n-1;i++) {
        paramRuler.push_back(i*1.0/(2*n-2));
    }
    if (knot) {
        free(knot);
    }
    knotTemp.clear();
    knot = (float*)malloc(((paramRuler.size()+1)/2+4)*sizeof(float));
    knot[0] = 0;
    knot[1] = 0;
    knot[2] = 0;
    knot[3] = 0;
    knotTemp.push_back(0);
    knotTemp.push_back(0);
    knotTemp.push_back(0);
    knotTemp.push_back(0);
    for (int i=4;i<=(paramRuler.size()+1)/2-1;i++) {
        knot[i] = (i-3)*1.0/((paramRuler.size()+1)/2-3);
        knotTemp.push_back((i-3)*1.0/((paramRuler.size()+1)/2-3));
    }
    knot[(paramRuler.size()+1)/2] = 1;
    knot[(paramRuler.size()+1)/2+1] = 1;
    knot[(paramRuler.size()+1)/2+2] = 1;
    knot[(paramRuler.size()+1)/2+3] = 1;
    knotTemp.push_back(1);
    knotTemp.push_back(1);
    knotTemp.push_back(1);
    knotTemp.push_back(1);
    gridLength = (paramRuler.size()+1)/2-1;

    while (flagErrorControl) {
        iterationCount++;
        lGrid = 1.0/((paramRuler.size()+1)/2-1);
        cout << "TPMS2STEP > SchwarzP Error Control Iteration: " << iterationCount << endl;
        int length = (paramRuler.size()+1)/2;// the number of control points is the same as the number of sample points
        int halfLength = (length+1)/2;
        sampleSchwarzP(paramRuler);
        pythonEXEC(this->modelType, this->getOffsetValue1(), this->getOffsetValue2());
        // Define the number of control points in the U and V directions.
        (*NumRows) = (paramRuler.size()+1)/2;
        (*NumCols) = (paramRuler.size()+1)/2;
        controlPointsPlusPIA.resize((*NumRows)*(*NumCols), 3);
        controlPointsMinusPIA.resize((*NumRows)*(*NumCols), 3);
        *numberOfPointCloud = (paramRuler.size()+1)/2*(paramRuler.size()+1)/2;
        U.resize(*numberOfPointCloud);
        V.resize(*numberOfPointCloud);
        surfacePlusMatrix.resize(*numberOfPointCloud, 3);
        surfaceMinusMatrix.resize(*numberOfPointCloud, 3);

        // read sample points
        readSemiFundamentalPatch(surfacePlusMatrix, surfaceMinusMatrix, numberOfPointCloud, this->getOffsetValue1());
        ofs.open("../data/uv/uv_SchwarzP_pia.off", ios::out);
        ofs << "OFF\n" << to_string(*numberOfPointCloud) << " 0 0" << endl;
        for (int i=0;i<(paramRuler.size()+1)/2;i++) {
            for (int j=0;j<(paramRuler.size()+1)/2;j++) {
                U(i*(paramRuler.size()+1)/2+j) = j*1.0/((paramRuler.size()+1)/2-1);
                V(i*(paramRuler.size()+1)/2+j) = i*1.0/((paramRuler.size()+1)/2-1);
                ofs << j*1.0/((paramRuler.size()+1)/2-1) << " " << i*1.0/((paramRuler.size()+1)/2-1) << " 0" << endl;
            }
        }
        ofs.close();
        cout << "TPMS2STEP > GPU constrained-PIA" << endl;
        constrainedPIAGPU(knot, 3, (*numberOfPointCloud), (*NumRows), (*NumCols), surfacePlusMatrix, surfaceMinusMatrix, U, V,
        &controlPointsPlusPIA, &controlPointsMinusPIA, this->getOffsetValue1(), this->modelType);

        // evaluate the error
        gridLength = length-1; // the grid number of control points in u or v direction
        maxNURBSerror = -1;
        std::vector<int> tmpRef;
        for (int i=0;i<gridLength;i++) {
            for (int j=0;j<gridLength;j++) {
                float aijx = controlPointsPlusPIA(i*length+j,0)+controlPointsPlusPIA((i+1)*length+(j+1),0)
                -controlPointsPlusPIA((i+1)*length+j,0)-controlPointsPlusPIA(i*length+(j+1),0);
                float aijy = controlPointsPlusPIA(i*length+j,1)+controlPointsPlusPIA((i+1)*length+(j+1),1)
                -controlPointsPlusPIA((i+1)*length+j,1)-controlPointsPlusPIA(i*length+(j+1),1);
                float aijz = controlPointsPlusPIA(i*length+j,2)+controlPointsPlusPIA((i+1)*length+(j+1),2)
                -controlPointsPlusPIA((i+1)*length+j,2)-controlPointsPlusPIA(i*length+(j+1),2);
                float error = (lGrid*lGrid*gridLength*gridLength*(std::sqrt(aijx*aijx+aijy*aijy+aijz*aijz)))/4;
                
                if (error > epsilon_2) {
                    // refine
                    tmpRef.push_back(j);
                    tmpRef.push_back(gridLength-1-j);
                    tmpRef.push_back(i);
                    tmpRef.push_back(gridLength-1-i);
                    sort(tmpRef.begin(), tmpRef.end());
                    tmpRef.erase(unique(tmpRef.begin(), tmpRef.end()), tmpRef.end());
                }
                if (maxNURBSerror < error) {
                    maxNURBSerror = error;
                }
            }
        }
        for (int i=0;i<tmpRef.size();i++) {
            refineTimes[tmpRef[i]] = 1;
        }
        if (maxNURBSerror <= epsilon_2 || iterationCount>=maxIteration) {
            flagErrorControl = false;
            cout << "TPMS2STEP > Error Control Iteration over." << endl;
            break;
        }
        int totalRefineNumber = paramRuler.size()+1;
        paramRuler.clear();
        for (int i=0;i<totalRefineNumber;i++) {
            paramRuler.push_back(i*1.0/(totalRefineNumber-1));
        }
        refineTimes.clear();
        refineTimes.resize(paramRuler.size()-1);
        for (int i=0;i<paramRuler.size()-1;i++) {
            refineTimes.push_back(0);
        }
        if (knot) {
            free(knot);
            knotTemp.clear();
            knot = (float*)malloc(((paramRuler.size()+1)/2+4)*sizeof(float));
            knot[0] = 0;
            knot[1] = 0;
            knot[2] = 0;
            knot[3] = 0;
            knotTemp.push_back(0);
            knotTemp.push_back(0);
            knotTemp.push_back(0);
            knotTemp.push_back(0);
            for (int i=4;i<=(paramRuler.size()+1)/2-1;i++) {
                knot[i] = (i-3)*1.0/((paramRuler.size()+1)/2-3);
                knotTemp.push_back((i-3)*1.0/((paramRuler.size()+1)/2-3));
            }
            knot[(paramRuler.size()+1)/2] = 1;
            knot[(paramRuler.size()+1)/2+1] = 1;
            knot[(paramRuler.size()+1)/2+2] = 1;
            knot[(paramRuler.size()+1)/2+3] = 1;
            knotTemp.push_back(1);
            knotTemp.push_back(1);
            knotTemp.push_back(1);
            knotTemp.push_back(1);
        }
    }
}

void TPMS_P::readSemiFundamentalPatch(Eigen::MatrixXf& surfacePlusMatrix, Eigen::MatrixXf& surfaceMinusMatrix, int* numberOfPointCloud, float offsetValue) {
    string tempString;
    ifstream ifs;
    int vertex, edge, face;
    float a, b, c;
    ifs.open("../data/SFP/SchwarzP/SFP_SchwarzP_offset_" + to_string(int(offsetValue * (-100) + 1 - 1)) + ".off", ios::in);
    ifs >> tempString;
    ifs >> vertex >> edge >> face;
    (*numberOfPointCloud) = vertex;
    for (Standard_Integer i = 0; i < vertex; i++) {
        ifs >> a >> b >> c;
        surfacePlusMatrix(i, 0) = a;
        surfacePlusMatrix(i, 1) = b;
        surfacePlusMatrix(i, 2) = c;
    }
    ifs.close();
    ifs.open("../data/SFP/SchwarzP/SFP_SchwarzP_offset_" + to_string(int((offsetValue * (-100) + 1 - 1) * (-1))) + ".off", ios::in);
    ifs >> tempString;
    ifs >> vertex >> edge >> face;
    for (Standard_Integer i = 0; i < vertex; i++) {
        ifs >> a >> b >> c;
        surfaceMinusMatrix(i, 0) = a;
        surfaceMinusMatrix(i, 1) = b;
        surfaceMinusMatrix(i, 2) = c;
    }
    ifs.close();
}

void TPMS_P::buildCurvesSchwarzP(std::vector<Handle(Geom_BSplineCurve)>& newCurve, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector, std::vector<Geom_BSplineSurface>& basicSurfaceElementVector) {
    newCurve.push_back(bSplineCurveVector[0]);
    newCurve.push_back(bSplineCurveVector[2]);

    newCurve.push_back(bSplineCurveVector[3]);
    newCurve.push_back(bSplineCurveVector[4]);

    newCurve.push_back(bSplineCurveVector[1]);
    newCurve.push_back(bSplineCurveVector[5]);

    newCurve.push_back(bSplineCurveVector[6]);
    newCurve.push_back(bSplineCurveVector[8]);

    newCurve.push_back(bSplineCurveVector[9]);
    newCurve.push_back(bSplineCurveVector[10]);

    newCurve.push_back(bSplineCurveVector[7]);
    newCurve.push_back(bSplineCurveVector[11]);

    std::vector<GeomFill_BSplineCurves> boundarySurfaceVector;
    GeomFill_BSplineCurves aFillSurface;

    Handle(Geom_BSplineCurve) copyCurveFirst = Handle(Geom_BSplineCurve)::DownCast(newCurve[0]->Copy());
    Handle(Geom_BSplineCurve) copyCurveSecond = Handle(Geom_BSplineCurve)::DownCast(newCurve[6]->Copy());
    gp_Pnt startPointFirst = copyCurveFirst->StartPoint();
    gp_Pnt startPointSecond = copyCurveSecond->StartPoint();
    gp_Pnt endPointFirst = copyCurveFirst->EndPoint();
    gp_Pnt endPointSecond = copyCurveSecond->EndPoint();

    TColgp_Array1OfPnt curveControl(1, 2);
    curveControl.SetValue(1, startPointFirst);
    curveControl.SetValue(2, startPointSecond);
    TColStd_Array1OfReal curveKnot(1, 2);
    curveKnot.SetValue(1, 0);
    curveKnot.SetValue(2, 1);
    TColStd_Array1OfInteger multiplicitiesCurve(1, 2);
    multiplicitiesCurve.SetValue(1, 2);
    multiplicitiesCurve.SetValue(2, 2);
    Handle(Geom_BSplineCurve) curveThird = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );

    curveControl.SetValue(1, endPointFirst);
    curveControl.SetValue(2, endPointSecond);
    Handle(Geom_BSplineCurve) curveFourth = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    aFillSurface = GeomFill_BSplineCurves(copyCurveFirst, copyCurveSecond, curveThird, curveFourth, GeomFill_CurvedStyle);
    boundarySurfaceVector.push_back(aFillSurface);

    copyCurveFirst = Handle(Geom_BSplineCurve)::DownCast(newCurve[1]->Copy());
    copyCurveSecond = Handle(Geom_BSplineCurve)::DownCast(newCurve[7]->Copy());

    startPointFirst = copyCurveFirst->StartPoint();
    startPointSecond = copyCurveSecond->StartPoint();
    endPointFirst = copyCurveFirst->EndPoint();
    endPointSecond = copyCurveSecond->EndPoint();

    curveControl.SetValue(1, startPointFirst);
    curveControl.SetValue(2, startPointSecond);
    curveThird = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    curveControl.SetValue(1, endPointFirst);
    curveControl.SetValue(2, endPointSecond);
    curveFourth = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    aFillSurface = GeomFill_BSplineCurves(copyCurveFirst, copyCurveSecond, curveThird, curveFourth, GeomFill_CurvedStyle);
    boundarySurfaceVector.push_back(aFillSurface);

    copyCurveFirst = Handle(Geom_BSplineCurve)::DownCast(newCurve[2]->Copy());
    copyCurveSecond = Handle(Geom_BSplineCurve)::DownCast(newCurve[8]->Copy());

    startPointFirst = copyCurveFirst->StartPoint();
    startPointSecond = copyCurveSecond->StartPoint();
    endPointFirst = copyCurveFirst->EndPoint();
    endPointSecond = copyCurveSecond->EndPoint();

    curveControl.SetValue(1, startPointFirst);
    curveControl.SetValue(2, startPointSecond);
    curveThird = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    curveControl.SetValue(1, endPointFirst);
    curveControl.SetValue(2, endPointSecond);
    curveFourth = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    aFillSurface = GeomFill_BSplineCurves(copyCurveFirst, copyCurveSecond, curveThird, curveFourth, GeomFill_CurvedStyle);
    boundarySurfaceVector.push_back(aFillSurface);

    copyCurveFirst = Handle(Geom_BSplineCurve)::DownCast(newCurve[3]->Copy());
    copyCurveSecond = Handle(Geom_BSplineCurve)::DownCast(newCurve[9]->Copy());

    startPointFirst = copyCurveFirst->StartPoint();
    startPointSecond = copyCurveSecond->StartPoint();
    endPointFirst = copyCurveFirst->EndPoint();
    endPointSecond = copyCurveSecond->EndPoint();

    curveControl.SetValue(1, startPointFirst);
    curveControl.SetValue(2, startPointSecond);
    curveThird = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    curveControl.SetValue(1, endPointFirst);
    curveControl.SetValue(2, endPointSecond);
    curveFourth = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    aFillSurface = GeomFill_BSplineCurves(copyCurveFirst, copyCurveSecond, curveThird, curveFourth, GeomFill_CurvedStyle);
    boundarySurfaceVector.push_back(aFillSurface);

    copyCurveFirst = Handle(Geom_BSplineCurve)::DownCast(newCurve[4]->Copy());
    copyCurveSecond = Handle(Geom_BSplineCurve)::DownCast(newCurve[10]->Copy());

    startPointFirst = copyCurveFirst->StartPoint();
    startPointSecond = copyCurveSecond->StartPoint();
    endPointFirst = copyCurveFirst->EndPoint();
    endPointSecond = copyCurveSecond->EndPoint();

    curveControl.SetValue(1, startPointFirst);
    curveControl.SetValue(2, startPointSecond);
    curveThird = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    curveControl.SetValue(1, endPointFirst);
    curveControl.SetValue(2, endPointSecond);
    curveFourth = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    aFillSurface = GeomFill_BSplineCurves(copyCurveFirst, copyCurveSecond, curveThird, curveFourth, GeomFill_CurvedStyle);
    boundarySurfaceVector.push_back(aFillSurface);

    copyCurveFirst = Handle(Geom_BSplineCurve)::DownCast(newCurve[5]->Copy());
    copyCurveSecond = Handle(Geom_BSplineCurve)::DownCast(newCurve[11]->Copy());

    startPointFirst = copyCurveFirst->StartPoint();
    startPointSecond = copyCurveSecond->StartPoint();
    endPointFirst = copyCurveFirst->EndPoint();
    endPointSecond = copyCurveSecond->EndPoint();

    curveControl.SetValue(1, startPointFirst);
    curveControl.SetValue(2, startPointSecond);
    curveThird = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    curveControl.SetValue(1, endPointFirst);
    curveControl.SetValue(2, endPointSecond);
    curveFourth = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    aFillSurface = GeomFill_BSplineCurves(copyCurveFirst, copyCurveSecond, curveThird, curveFourth, GeomFill_CurvedStyle);
    boundarySurfaceVector.push_back(aFillSurface);

    gp_Ax2 xA(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)); // planar symetry
    gp_Ax2 yA(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)); // planar symetry
    gp_Ax2 zA(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)); // planar symetry
    gp_Trsf xFace, yFace, zFace;
    xFace.SetMirror(xA);
    yFace.SetMirror(yA);
    zFace.SetMirror(zA);

    Handle(Geom_BSplineSurface) newSurface1 = boundarySurfaceVector[0].Surface();
    newSurface1->Translate(gp_Vec(0, 0, (-1)* cellSizeZ));
    Handle(Geom_BSplineSurface) newSurface2 = boundarySurfaceVector[1].Surface();
    newSurface2->Translate(gp_Vec(0, 0, (-1)* cellSizeZ));

    Handle(Geom_BSplineSurface) newSurface3 = Handle(Geom_BSplineSurface)::DownCast(newSurface1->Transformed(xFace));
    Handle(Geom_BSplineSurface) newSurface4 = Handle(Geom_BSplineSurface)::DownCast(newSurface2->Transformed(xFace));
    Handle(Geom_BSplineSurface) newSurface5 = Handle(Geom_BSplineSurface)::DownCast(newSurface1->Transformed(yFace));
    Handle(Geom_BSplineSurface) newSurface6 = Handle(Geom_BSplineSurface)::DownCast(newSurface2->Transformed(yFace));
    Handle(Geom_BSplineSurface) newSurface7 = Handle(Geom_BSplineSurface)::DownCast(newSurface5->Transformed(xFace));
    Handle(Geom_BSplineSurface) newSurface8 = Handle(Geom_BSplineSurface)::DownCast(newSurface6->Transformed(xFace));

    basicSurfaceElementVector.push_back(*newSurface1);
    basicSurfaceElementVector.push_back(*newSurface2);
    basicSurfaceElementVector.push_back(*newSurface3);
    basicSurfaceElementVector.push_back(*newSurface4);
    basicSurfaceElementVector.push_back(*newSurface5);
    basicSurfaceElementVector.push_back(*newSurface6);
    basicSurfaceElementVector.push_back(*newSurface7);
    basicSurfaceElementVector.push_back(*newSurface8);

    newSurface1 = boundarySurfaceVector[2].Surface();
    newSurface1->Translate(gp_Vec((-1)* cellSizeX, 0, 0));
    newSurface2 = boundarySurfaceVector[3].Surface();
    newSurface2->Translate(gp_Vec((-1)* cellSizeX, 0, 0));

    newSurface3 = Handle(Geom_BSplineSurface)::DownCast(newSurface1->Transformed(yFace));
    newSurface4 = Handle(Geom_BSplineSurface)::DownCast(newSurface2->Transformed(yFace));
    newSurface5 = Handle(Geom_BSplineSurface)::DownCast(newSurface1->Transformed(zFace));
    newSurface6 = Handle(Geom_BSplineSurface)::DownCast(newSurface2->Transformed(zFace));
    newSurface7 = Handle(Geom_BSplineSurface)::DownCast(newSurface5->Transformed(yFace));
    newSurface8 = Handle(Geom_BSplineSurface)::DownCast(newSurface6->Transformed(yFace));

    basicSurfaceElementVector.push_back(*newSurface1);
    basicSurfaceElementVector.push_back(*newSurface2);
    basicSurfaceElementVector.push_back(*newSurface3);
    basicSurfaceElementVector.push_back(*newSurface4);
    basicSurfaceElementVector.push_back(*newSurface5);
    basicSurfaceElementVector.push_back(*newSurface6);
    basicSurfaceElementVector.push_back(*newSurface7);
    basicSurfaceElementVector.push_back(*newSurface8);

    newSurface1 = boundarySurfaceVector[4].Surface();
    newSurface1->Translate(gp_Vec(0, (-1)* cellSizeY, 0));
    newSurface2 = boundarySurfaceVector[5].Surface();
    newSurface2->Translate(gp_Vec(0, (-1)* cellSizeY, 0));

    newSurface3 = Handle(Geom_BSplineSurface)::DownCast(newSurface1->Transformed(xFace));
    newSurface4 = Handle(Geom_BSplineSurface)::DownCast(newSurface2->Transformed(xFace));
    newSurface5 = Handle(Geom_BSplineSurface)::DownCast(newSurface1->Transformed(zFace));
    newSurface6 = Handle(Geom_BSplineSurface)::DownCast(newSurface2->Transformed(zFace));
    newSurface7 = Handle(Geom_BSplineSurface)::DownCast(newSurface5->Transformed(xFace));
    newSurface8 = Handle(Geom_BSplineSurface)::DownCast(newSurface6->Transformed(xFace));

    basicSurfaceElementVector.push_back(*newSurface1);
    basicSurfaceElementVector.push_back(*newSurface2);
    basicSurfaceElementVector.push_back(*newSurface3);
    basicSurfaceElementVector.push_back(*newSurface4);
    basicSurfaceElementVector.push_back(*newSurface5);
    basicSurfaceElementVector.push_back(*newSurface6);
    basicSurfaceElementVector.push_back(*newSurface7);
    basicSurfaceElementVector.push_back(*newSurface8);
}