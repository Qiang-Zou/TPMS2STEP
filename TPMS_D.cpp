#include "TPMS_D.h"

TPMS_D::TPMS_D()
{
}

TPMS_D::~TPMS_D()
{
}

bool TPMS_D::constructSolidModel() {
    this->cellScaleX = 1.0;
    this->cellScaleY = 1.0;
    this->cellScaleZ = 1.0;
    this->cellSizeX = 4.0 * this->cellScaleX;
    this->cellSizeY = 4.0 * this->cellScaleY;
    this->cellSizeZ = 4.0 * this->cellScaleZ;
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
    int numberOfQuadMesh = 6;
    std::vector<float> quadMeshPoints_1;
    std::vector<float> quadMeshPoints_2;
    std::vector<float> quadMeshPointsTotal;
#endif

    std::vector<math_Matrix> surfaceMatrixVector;
    std::vector<math_Matrix> surfaceMatrixVectorSecond;
    std::vector<Handle(Geom_BSplineSurface)> cellSurfaceVector;
    std::vector<TColgp_Array2OfPnt> cellSurfacePointsVector;
    for (int i=1;i<=NumRows;i++) {
        for (int j=1;j<=NumCols;j++) {
            (*ctrlpts_1).SetValue(i,j,gp_Pnt(controlPointsPlusPIA((i-1)*NumCols+(j-1), 0), controlPointsPlusPIA((i-1)*NumCols+(j-1), 1), controlPointsPlusPIA((i-1)*NumCols+(j-1), 2)));
        }
    }
    for (int i=1;i<=NumRows;i++) {
        for (int j=1;j<=NumCols;j++) {
            (*ctrlpts_2).SetValue(i,j,gp_Pnt(controlPointsMinusPIA((i-1)*NumCols+(j-1), 0), controlPointsMinusPIA((i-1)*NumCols+(j-1), 1), controlPointsMinusPIA((i-1)*NumCols+(j-1), 2)));
        }
    }

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
#endif

    int controlPointsColLength = (*ctrlpts_1).ColLength();
    int controlPointsRowLength = (*ctrlpts_1).RowLength();

    math_Matrix surface_0(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_1(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_2(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_3(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_4(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_5(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_6(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_0_Second(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_1_Second(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_2_Second(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_3_Second(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_4_Second(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_5_Second(1, (*ctrlpts_1).Size(), 1, 3);
    math_Matrix surface_6_Second(1, (*ctrlpts_1).Size(), 1, 3);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadMeshMatrix_1(1, numberOfQuadMesh*numberOfQuadMesh, 1, 3);
    math_Matrix quadMeshMatrix_2(1, numberOfQuadMesh*numberOfQuadMesh, 1, 3);

    std::vector<math_Matrix> quadMeshVector_1;
    std::vector<math_Matrix> quadMeshVector_2;

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

    TColgp_Array2OfPnt controlPoints_1(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_2(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_3(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_4(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_5(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_6(1, controlPointsColLength, 1, controlPointsRowLength);

    double ksi = 1.0 / sqrt(2.0);
    math_Matrix F1(1, 3, 1, 3, 0);
    math_Matrix F2(1, 3, 1, 3, 0);
    math_Matrix F3(1, 3, 1, 3, 0);
    math_Matrix F4(1, 3, 1, 3, 0);
    math_Matrix F5(1, 3, 1, 3, 0);
    math_Matrix F6(1, 3, 1, 3, 0);
    F1(1, 1) = ksi; F1(1, 2) = (-1) * ksi; F1(1, 3) = 0;
    F1(2, 1) = ksi; F1(2, 2) = ksi; F1(2, 3) = 0;
    F1(3, 1) = 0; F1(3, 2) = 0; F1(3, 3) = 1;

    F2(1, 1) = (-1) * ksi; F2(1, 2) = ksi; F2(1, 3) = 0;
    F2(2, 1) = 0; F2(2, 2) = 0; F2(2, 3) = 1;
    F2(3, 1) = ksi; F2(3, 2) = ksi; F2(3, 3) = 0;

    F3(1, 1) = 0; F3(1, 2) = 0; F3(1, 3) = 1;
    F3(2, 1) = (-1) * ksi; F3(2, 2) = (-1) * ksi; F3(2, 3) = 0;
    F3(3, 1) = ksi; F3(3, 2) = (-1) * ksi; F3(3, 3) = 0;

    F4(1, 1) = (-1) * ksi; F4(1, 2) = (-1) * ksi; F4(1, 3) = 0;
    F4(2, 1) = (-1) * ksi; F4(2, 2) = ksi; F4(2, 3) = 0;
    F4(3, 1) = 0; F4(3, 2) = 0; F4(3, 3) = -1;

    F5(1, 1) = 0; F5(1, 2) = 0; F5(1, 3) = -1;
    F5(2, 1) = ksi; F5(2, 2) = (-1) * ksi; F5(2, 3) = 0;
    F5(3, 1) = (-1) * ksi; F5(3, 2) = (-1) * ksi; F5(3, 3) = 0;

    F6(1, 1) = ksi; F6(1, 2) = ksi; F6(1, 3) = 0;
    F6(2, 1) = 0; F6(2, 2) = 0; F6(2, 3) = -1;
    F6(3, 1) = (-1) * ksi; F6(3, 2) = ksi; F6(3, 3) = 0;

    std::vector<math_Matrix> Q;
    std::vector<math_Vector> eta;
    math_Matrix Q1(1, 3, 1, 3, 0);
    math_Matrix Q2(1, 3, 1, 3, 0);
    math_Matrix Q3(1, 3, 1, 3, 0);
    math_Matrix Q4(1, 3, 1, 3, 0);
    math_Vector eta1(1, 3, 0);
    math_Vector eta2(1, 3, 0);
    math_Vector eta3(1, 3, 0);
    math_Vector eta4(1, 3, 0);

    Q1(1, 1) = 0; Q1(1, 2) = -1; Q1(1, 3) = 0;
    Q1(2, 1) = -1; Q1(2, 2) = 0; Q1(2, 3) = 0;
    Q1(3, 1) = 0; Q1(3, 2) = 0; Q1(3, 3) = 1;
    eta1(1) = 1; eta1(2) = 1; eta1(3) = 0;

    Q2(1, 1) = 0; Q2(1, 2) = 1; Q2(1, 3) = 0;
    Q2(2, 1) = -1; Q2(2, 2) = 0; Q2(2, 3) = 0;
    Q2(3, 1) = 0; Q2(3, 2) = 0; Q2(3, 3) = 1;
    eta2(1) = 0; eta2(2) = 0; eta2(3) = -1;

    Q3(1, 1) = 0; Q3(1, 2) = 1; Q3(1, 3) = 0;
    Q3(2, 1) = 1; Q3(2, 2) = 0; Q3(2, 3) = 0;
    Q3(3, 1) = 0; Q3(3, 2) = 0; Q3(3, 3) = 1;
    eta3(1) = -1; eta3(2) = -1; eta3(3) = 0;

    Q4(1, 1) = 0; Q4(1, 2) = -1; Q4(1, 3) = 0;
    Q4(2, 1) = 1; Q4(2, 2) = 0; Q4(2, 3) = 0;
    Q4(3, 1) = 0; Q4(3, 2) = 0; Q4(3, 3) = 1;
    eta4(1) = 0; eta4(2) = 0; eta4(3) = -1;

    Q.push_back(Q1);
    Q.push_back(Q2);
    Q.push_back(Q3);
    Q.push_back(Q4);

    eta.push_back(eta1);
    eta.push_back(eta2);
    eta.push_back(eta3);
    eta.push_back(eta4);

    for (int i = 1; i <= (*ctrlpts_1).ColLength(); i++) {
        for (int j = 1; j <= (*ctrlpts_1).RowLength(); j++) {
            surface_0((i - 1) * controlPointsRowLength + j, 1) = (*ctrlpts_1).Value(i, j).X();
            surface_0((i - 1) * controlPointsRowLength + j, 2) = (*ctrlpts_1).Value(i, j).Y();
            surface_0((i - 1) * controlPointsRowLength + j, 3) = (*ctrlpts_1).Value(i, j).Z();
            // cout << (*ctrlpts_1).Value(i, j).X() << " " << (*ctrlpts_1).Value(i, j).Y() << " " << (*ctrlpts_1).Value(i, j).Z() << endl;
        }
    }
    for (int i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (int j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_0_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1) = (*ctrlpts_2).Value(i, j).X();
            surface_0_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2) = (*ctrlpts_2).Value(i, j).Y();
            surface_0_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) = (*ctrlpts_2).Value(i, j).Z();
            // cout << (*ctrlpts_2).Value(i, j).X() << " " << (*ctrlpts_2).Value(i, j).Y() << " " << (*ctrlpts_2).Value(i, j).Z() << endl;
        }
    }

    math_Matrix temp1 = F1 * surface_0.Transposed();
    surface_1 = temp1.Transposed();
    surfaceMatrixVector.push_back(surface_1);
    for (Standard_Integer i = 1; i <= (*ctrlpts_1).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_1).RowLength(); j++) {
            controlPoints_1.SetValue(i, j, gp_Pnt(surface_1((i - 1) * controlPointsRowLength + j, 1), surface_1((i - 1) * controlPointsRowLength + j, 2), surface_1((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_1);
    
#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp1 = F1 * quadMeshMatrix_1.Transposed();
    math_Matrix quadTemp12 = quadTemp1.Transposed();
    quadMeshVector_1.push_back(quadTemp12);
#endif

    math_Matrix temp2 = F2 * surface_0_Second.Transposed();
    math_Matrix temp22 = temp2.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_2((i - 1) * controlPointsRowLength + j, 1) = temp22((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_2((i - 1) * controlPointsRowLength + j, 2) = temp22((i - 1) * controlPointsRowLength + j, 2);
            surface_2((i - 1) * controlPointsRowLength + j, 3) = temp22((i - 1) * controlPointsRowLength + j, 3);
        }
    }
    surfaceMatrixVector.push_back(surface_2);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_2.SetValue(i, j, gp_Pnt(surface_2((i - 1) * controlPointsRowLength + j, 1), surface_2((i - 1) * controlPointsRowLength + j, 2), surface_2((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_2);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp2 = F2 * quadMeshMatrix_2.Transposed();
    math_Matrix quadTemp22 = quadTemp2.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp22(i*numberOfQuadMesh+j+1, 1) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp22);
#endif

    math_Matrix temp3 = F3 * surface_0_Second.Transposed();
    math_Matrix temp33 = temp3.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_3((i - 1) * controlPointsRowLength + j, 1) = temp33((i - 1) * controlPointsRowLength + j, 1);
            surface_3((i - 1) * controlPointsRowLength + j, 2) = temp33((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_3((i - 1) * controlPointsRowLength + j, 3) = temp33((i - 1) * controlPointsRowLength + j, 3);
        }
    }
    surfaceMatrixVector.push_back(surface_3);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_3.SetValue(i, j, gp_Pnt(surface_3((i - 1) * controlPointsRowLength + j, 1), surface_3((i - 1) * controlPointsRowLength + j, 2), surface_3((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_3);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp3 = F3 * quadMeshMatrix_2.Transposed();
    math_Matrix quadTemp32 = quadTemp3.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp32(i*numberOfQuadMesh+j+1, 2) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp32);
#endif

    math_Matrix temp4 = F4 * surface_0_Second.Transposed();
    math_Matrix temp44 = temp4.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_4((i - 1) * controlPointsRowLength + j, 1) = temp44((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_4((i - 1) * controlPointsRowLength + j, 2) = temp44((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_4((i - 1) * controlPointsRowLength + j, 3) = temp44((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_4);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_4.SetValue(i, j, gp_Pnt(surface_4((i - 1) * controlPointsRowLength + j, 1), surface_4((i - 1) * controlPointsRowLength + j, 2), surface_4((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_4);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp4 = F4 * quadMeshMatrix_2.Transposed();
    math_Matrix quadTemp42 = quadTemp4.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp42(i*numberOfQuadMesh+j+1, 1) += 1;
            quadTemp42(i*numberOfQuadMesh+j+1, 2) += 1;
            quadTemp42(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp42);
#endif

    math_Matrix temp5 = F5 * surface_0.Transposed();
    math_Matrix temp55 = temp5.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_5((i - 1) * controlPointsRowLength + j, 1) = temp55((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_5((i - 1) * controlPointsRowLength + j, 2) = temp55((i - 1) * controlPointsRowLength + j, 2);
            surface_5((i - 1) * controlPointsRowLength + j, 3) = temp55((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_5);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_5.SetValue(i, j, gp_Pnt(surface_5((i - 1) * controlPointsRowLength + j, 1), surface_5((i - 1) * controlPointsRowLength + j, 2), surface_5((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_5);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp5 = F5 * quadMeshMatrix_1.Transposed();
    math_Matrix quadTemp52 = quadTemp5.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp52(i*numberOfQuadMesh+j+1, 1) += 1;
            quadTemp52(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp52);
#endif

    math_Matrix temp6 = F6 * surface_0.Transposed();
    math_Matrix temp66 = temp6.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_6((i - 1) * controlPointsRowLength + j, 1) = temp66((i - 1) * controlPointsRowLength + j, 1);
            surface_6((i - 1) * controlPointsRowLength + j, 2) = temp66((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_6((i - 1) * controlPointsRowLength + j, 3) = temp66((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_6);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_6.SetValue(i, j, gp_Pnt(surface_6((i - 1) * controlPointsRowLength + j, 1), surface_6((i - 1) * controlPointsRowLength + j, 2), surface_6((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_6);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp6 = F6 * quadMeshMatrix_1.Transposed();
    math_Matrix quadTemp62 = quadTemp6.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp62(i*numberOfQuadMesh+j+1, 2) += 1;
            quadTemp62(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp62);
#endif

    temp1 = F1 * surface_0_Second.Transposed();
    surface_1_Second = temp1.Transposed();
    surfaceMatrixVectorSecond.push_back(surface_1_Second);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_1.SetValue(i, j, gp_Pnt(surface_1_Second((i - 1) * controlPointsRowLength + j, 1), surface_1_Second((i - 1) * controlPointsRowLength + j, 2), surface_1_Second((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_1);

#ifdef USE_QUAD_HEX_MESH
    quadTemp1 = F1 * quadMeshMatrix_2.Transposed();
    quadTemp12 = quadTemp1.Transposed();
    quadMeshVector_2.push_back(quadTemp12);
#endif

    temp2 = F2 * surface_0.Transposed();
    temp22 = temp2.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_2_Second((i - 1) * controlPointsRowLength + j, 1) = temp22((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_2_Second((i - 1) * controlPointsRowLength + j, 2) = temp22((i - 1) * controlPointsRowLength + j, 2);
            surface_2_Second((i - 1) * controlPointsRowLength + j, 3) = temp22((i - 1) * controlPointsRowLength + j, 3);
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_2_Second);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_2.SetValue(i, j, gp_Pnt(surface_2_Second((i - 1) * controlPointsRowLength + j, 1), surface_2_Second((i - 1) * controlPointsRowLength + j, 2), surface_2_Second((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_2);

#ifdef USE_QUAD_HEX_MESH
    quadTemp2 = F2 * quadMeshMatrix_1.Transposed();
    quadTemp22 = quadTemp2.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp22(i*numberOfQuadMesh+j+1, 1) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp22);
#endif

    temp3 = F3 * surface_0.Transposed();
    temp33 = temp3.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_3_Second((i - 1) * controlPointsRowLength + j, 1) = temp33((i - 1) * controlPointsRowLength + j, 1);
            surface_3_Second((i - 1) * controlPointsRowLength + j, 2) = temp33((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_3_Second((i - 1) * controlPointsRowLength + j, 3) = temp33((i - 1) * controlPointsRowLength + j, 3);
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_3_Second);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_3.SetValue(i, j, gp_Pnt(surface_3_Second((i - 1) * controlPointsRowLength + j, 1), surface_3_Second((i - 1) * controlPointsRowLength + j, 2), surface_3_Second((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_3);

#ifdef USE_QUAD_HEX_MESH
    quadTemp3 = F3 * quadMeshMatrix_1.Transposed();
    quadTemp32 = quadTemp3.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp32(i*numberOfQuadMesh+j+1, 2) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp32);
#endif

    temp4 = F4 * surface_0.Transposed();
    temp44 = temp4.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_4_Second((i - 1) * controlPointsRowLength + j, 1) = temp44((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_4_Second((i - 1) * controlPointsRowLength + j, 2) = temp44((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_4_Second((i - 1) * controlPointsRowLength + j, 3) = temp44((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_4_Second);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_4.SetValue(i, j, gp_Pnt(surface_4_Second((i - 1) * controlPointsRowLength + j, 1), surface_4_Second((i - 1) * controlPointsRowLength + j, 2), surface_4_Second((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_4);

#ifdef USE_QUAD_HEX_MESH
    quadTemp4 = F4 * quadMeshMatrix_1.Transposed();
    quadTemp42 = quadTemp4.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp42(i*numberOfQuadMesh+j+1, 1) += 1;
            quadTemp42(i*numberOfQuadMesh+j+1, 2) += 1;
            quadTemp42(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp42);
#endif

    temp5 = F5 * surface_0_Second.Transposed();
    temp55 = temp5.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_5_Second((i - 1) * controlPointsRowLength + j, 1) = temp55((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_5_Second((i - 1) * controlPointsRowLength + j, 2) = temp55((i - 1) * controlPointsRowLength + j, 2);
            surface_5_Second((i - 1) * controlPointsRowLength + j, 3) = temp55((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_5_Second);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_5.SetValue(i, j, gp_Pnt(surface_5_Second((i - 1) * controlPointsRowLength + j, 1), surface_5_Second((i - 1) * controlPointsRowLength + j, 2), surface_5_Second((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_5);

#ifdef USE_QUAD_HEX_MESH
    quadTemp5 = F5 * quadMeshMatrix_2.Transposed();
    quadTemp52 = quadTemp5.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp52(i*numberOfQuadMesh+j+1, 1) += 1;
            quadTemp52(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp52);
#endif

    temp6 = F6 * surface_0_Second.Transposed();
    temp66 = temp6.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_6_Second((i - 1) * controlPointsRowLength + j, 1) = temp66((i - 1) * controlPointsRowLength + j, 1);
            surface_6_Second((i - 1) * controlPointsRowLength + j, 2) = temp66((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_6_Second((i - 1) * controlPointsRowLength + j, 3) = temp66((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_6_Second);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            controlPoints_6.SetValue(i, j, gp_Pnt(surface_6_Second((i - 1) * controlPointsRowLength + j, 1), surface_6_Second((i - 1) * controlPointsRowLength + j, 2), surface_6_Second((i - 1) * controlPointsRowLength + j, 3)));
        }
    }
    cellSurfacePointsVector.push_back(controlPoints_6);

#ifdef USE_QUAD_HEX_MESH
    quadTemp6 = F6 * quadMeshMatrix_2.Transposed();
    quadTemp62 = quadTemp6.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp62(i*numberOfQuadMesh+j+1, 2) += 1;
            quadTemp62(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp62);
#endif

    // duplicatet the faces for 8 times
    std::vector<TColgp_Array2OfPnt> controlPointsTotalVector;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 0 || j == 2) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
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
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 0 || j == 2) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 1 || j == 3) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3);
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 1 || j == 3) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += (eta[j].Value(1)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 1 || j == 3) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3);
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 1 || j == 3) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += (eta[j].Value(2)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 0 || j == 2) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3);
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 0 || j == 2) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += (eta[j].Value(1)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += (eta[j].Value(2)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 1 || j == 3) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3) + 2;
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 1 || j == 3) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += (eta[j].Value(3)+2);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 0 || j == 2) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3) + 2;
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 0 || j == 2) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += (eta[j].Value(1)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += (eta[j].Value(3)+2);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 1 || j == 3) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3) + 2;
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 1 || j == 3) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += (eta[j].Value(1)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += (eta[j].Value(2)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += (eta[j].Value(3)+2);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 0 || j == 2) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3) + 2;
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 0 || j == 2) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += (eta[j].Value(2)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += (eta[j].Value(3)+2);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    std::vector<TColgp_Array1OfPnt> curveControlPointsVector;
    std::vector<Handle(Geom_BSplineCurve)> bSplineCurveVector;
    buildCurvesDiamond(Uknot, multiplicitiesU, controlPointsColLength, cellSurfacePointsVector, curveControlPointsVector, bSplineCurveVector, controlPointsTotalVector);
    
    surfaceMatrixVector.clear();
    surfaceMatrixVectorSecond.clear();

    // read sample points
    readSemiFundamentalPatch(surfacePlusMatrix, surfaceMinusMatrix, &numberOfPointCloud, this->getOffsetValue2());
    constrainedPIAGPU(knot, 3, numberOfPointCloud, NumRows, NumCols, surfacePlusMatrix, surfaceMinusMatrix, U, V,
    &controlPointsPlusPIA, &controlPointsMinusPIA, this->getOffsetValue2(), this->modelType);
    
    // ofs.open("controlPoint2.txt", ios::out);
    // for (int i=0;i<NumRows*NumCols;i++) {
    //     ofs << controlPointsPlusPIA(i,0) << " " << controlPointsPlusPIA(i,1) << " " << controlPointsPlusPIA(i,2) << endl;
    // }
    // ofs.close();
    
    for (int i=1;i<=NumRows;i++) {
        for (int j=1;j<=NumCols;j++) {
            (*ctrlpts_1).SetValue(i,j,gp_Pnt(controlPointsPlusPIA((i-1)*NumCols+(j-1), 0), controlPointsPlusPIA((i-1)*NumCols+(j-1), 1), controlPointsPlusPIA((i-1)*NumCols+(j-1), 2)));
        }
    }
    for (int i=1;i<=NumRows;i++) {
        for (int j=1;j<=NumCols;j++) {
            (*ctrlpts_2).SetValue(i,j,gp_Pnt(controlPointsMinusPIA((i-1)*NumCols+(j-1), 0), controlPointsMinusPIA((i-1)*NumCols+(j-1), 1), controlPointsMinusPIA((i-1)*NumCols+(j-1), 2)));
        }
    }

#ifdef USE_QUAD_HEX_MESH
    Handle(Geom_BSplineSurface) thirdSplineSurface = new Geom_BSplineSurface(
        *ctrlpts_1,
        *Uknot,
        *Vknot,
        *multiplicitiesU,
        *multiplicitiesV,
        3,
        3
    );
    Handle(Geom_BSplineSurface) fourthSplineSurface = new Geom_BSplineSurface(
        *ctrlpts_2,
        *Uknot,
        *Vknot,
        *multiplicitiesU,
        *multiplicitiesV,
        3,
        3
    );
    quadMeshPoints_1.clear();
    quadMeshPoints_2.clear();
    quadMeshVector_1.clear();
    quadMeshVector_2.clear();
    extractMesh(thirdSplineSurface, quadMeshPoints_1, numberOfQuadMesh);
    extractMesh(fourthSplineSurface, quadMeshPoints_2, numberOfQuadMesh);
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

    for (int i = 1; i <= controlPointsColLength; i++) {
        for (int j = 1; j <= controlPointsRowLength; j++) {
            surface_0((i - 1) * controlPointsRowLength + j, 1) = (*ctrlpts_1).Value(i, j).X();
            surface_0((i - 1) * controlPointsRowLength + j, 2) = (*ctrlpts_1).Value(i, j).Y();
            surface_0((i - 1) * controlPointsRowLength + j, 3) = (*ctrlpts_1).Value(i, j).Z();
        }
    }
    for (int i = 1; i <= controlPointsColLength; i++) {
        for (int j = 1; j <= controlPointsRowLength; j++) {
            surface_0_Second((i - 1) * controlPointsRowLength + j, 1) = (*ctrlpts_2).Value(i, j).X();
            surface_0_Second((i - 1) * controlPointsRowLength + j, 2) = (*ctrlpts_2).Value(i, j).Y();
            surface_0_Second((i - 1) * controlPointsRowLength + j, 3) = (*ctrlpts_2).Value(i, j).Z();
        }
    }

    temp1 = F1 * surface_0.Transposed();
    surface_1 = temp1.Transposed();
    surfaceMatrixVector.push_back(surface_1);

#ifdef USE_QUAD_HEX_MESH
    quadTemp1 = F1 * quadMeshMatrix_1.Transposed();
    quadTemp12 = quadTemp1.Transposed();
    quadMeshVector_1.push_back(quadTemp12);
#endif

    temp2 = F2 * surface_0_Second.Transposed();
    temp22 = temp2.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_2((i - 1) * controlPointsRowLength + j, 1) = temp22((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_2((i - 1) * controlPointsRowLength + j, 2) = temp22((i - 1) * controlPointsRowLength + j, 2);
            surface_2((i - 1) * controlPointsRowLength + j, 3) = temp22((i - 1) * controlPointsRowLength + j, 3);
        }
    }
    surfaceMatrixVector.push_back(surface_2);

#ifdef USE_QUAD_HEX_MESH
    quadTemp2 = F2 * quadMeshMatrix_2.Transposed();
    quadTemp22 = quadTemp2.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp22(i*numberOfQuadMesh+j+1, 1) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp22);
#endif

    temp3 = F3 * surface_0_Second.Transposed();
    temp33 = temp3.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_3((i - 1) * controlPointsRowLength + j, 1) = temp33((i - 1) * controlPointsRowLength + j, 1);
            surface_3((i - 1) * controlPointsRowLength + j, 2) = temp33((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_3((i - 1) * controlPointsRowLength + j, 3) = temp33((i - 1) * controlPointsRowLength + j, 3);
        }
    }
    surfaceMatrixVector.push_back(surface_3);

#ifdef USE_QUAD_HEX_MESH
    quadTemp3 = F3 * quadMeshMatrix_2.Transposed();
    quadTemp32 = quadTemp3.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp32(i*numberOfQuadMesh+j+1, 2) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp32);
#endif

    temp4 = F4 * surface_0_Second.Transposed();
    temp44 = temp4.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_4((i - 1) * controlPointsRowLength + j, 1) = temp44((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_4((i - 1) * controlPointsRowLength + j, 2) = temp44((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_4((i - 1) * controlPointsRowLength + j, 3) = temp44((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_4);

#ifdef USE_QUAD_HEX_MESH
    quadTemp4 = F4 * quadMeshMatrix_2.Transposed();
    quadTemp42 = quadTemp4.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp42(i*numberOfQuadMesh+j+1, 1) += 1;
            quadTemp42(i*numberOfQuadMesh+j+1, 2) += 1;
            quadTemp42(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp42);
#endif

    temp5 = F5 * surface_0.Transposed();
    temp55 = temp5.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_5((i - 1) * controlPointsRowLength + j, 1) = temp55((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_5((i - 1) * controlPointsRowLength + j, 2) = temp55((i - 1) * controlPointsRowLength + j, 2);
            surface_5((i - 1) * controlPointsRowLength + j, 3) = temp55((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_5);

#ifdef USE_QUAD_HEX_MESH
    quadTemp5 = F5 * quadMeshMatrix_1.Transposed();
    quadTemp52 = quadTemp5.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp52(i*numberOfQuadMesh+j+1, 1) += 1;
            quadTemp52(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp52);
#endif

    temp6 = F6 * surface_0.Transposed();
    temp66 = temp6.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_6((i - 1) * controlPointsRowLength + j, 1) = temp66((i - 1) * controlPointsRowLength + j, 1);
            surface_6((i - 1) * controlPointsRowLength + j, 2) = temp66((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_6((i - 1) * controlPointsRowLength + j, 3) = temp66((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_6);

#ifdef USE_QUAD_HEX_MESH
    quadTemp6 = F6 * quadMeshMatrix_1.Transposed();
    quadTemp62 = quadTemp6.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp62(i*numberOfQuadMesh+j+1, 2) += 1;
            quadTemp62(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp62);
#endif

    temp1 = F1 * surface_0_Second.Transposed();
    surface_1_Second = temp1.Transposed();
    surfaceMatrixVectorSecond.push_back(surface_1_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp1 = F1 * quadMeshMatrix_2.Transposed();
    quadTemp12 = quadTemp1.Transposed();
    quadMeshVector_2.push_back(quadTemp12);
#endif

    temp2 = F2 * surface_0.Transposed();
    temp22 = temp2.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_2_Second((i - 1) * controlPointsRowLength + j, 1) = temp22((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_2_Second((i - 1) * controlPointsRowLength + j, 2) = temp22((i - 1) * controlPointsRowLength + j, 2);
            surface_2_Second((i - 1) * controlPointsRowLength + j, 3) = temp22((i - 1) * controlPointsRowLength + j, 3);
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_2_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp2 = F2 * quadMeshMatrix_1.Transposed();
    quadTemp22 = quadTemp2.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp22(i*numberOfQuadMesh+j+1, 1) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp22);
#endif

    temp3 = F3 * surface_0.Transposed();
    temp33 = temp3.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_3_Second((i - 1) * controlPointsRowLength + j, 1) = temp33((i - 1) * controlPointsRowLength + j, 1);
            surface_3_Second((i - 1) * controlPointsRowLength + j, 2) = temp33((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_3_Second((i - 1) * controlPointsRowLength + j, 3) = temp33((i - 1) * controlPointsRowLength + j, 3);
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_3_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp3 = F3 * quadMeshMatrix_1.Transposed();
    quadTemp32 = quadTemp3.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp32(i*numberOfQuadMesh+j+1, 2) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp32);
#endif

    temp4 = F4 * surface_0.Transposed();
    temp44 = temp4.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_4_Second((i - 1) * controlPointsRowLength + j, 1) = temp44((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_4_Second((i - 1) * controlPointsRowLength + j, 2) = temp44((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_4_Second((i - 1) * controlPointsRowLength + j, 3) = temp44((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_4_Second);
    
#ifdef USE_QUAD_HEX_MESH
    quadTemp4 = F4 * quadMeshMatrix_1.Transposed();
    quadTemp42 = quadTemp4.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp42(i*numberOfQuadMesh+j+1, 1) += 1;
            quadTemp42(i*numberOfQuadMesh+j+1, 2) += 1;
            quadTemp42(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp42);
#endif

    temp5 = F5 * surface_0_Second.Transposed();
    temp55 = temp5.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_5_Second((i - 1) * controlPointsRowLength + j, 1) = temp55((i - 1) * controlPointsRowLength + j, 1) + 1;
            surface_5_Second((i - 1) * controlPointsRowLength + j, 2) = temp55((i - 1) * controlPointsRowLength + j, 2);
            surface_5_Second((i - 1) * controlPointsRowLength + j, 3) = temp55((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_5_Second);
    
#ifdef USE_QUAD_HEX_MESH
    quadTemp5 = F5 * quadMeshMatrix_2.Transposed();
    quadTemp52 = quadTemp5.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp52(i*numberOfQuadMesh+j+1, 1) += 1;
            quadTemp52(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp52);
#endif

    temp6 = F6 * surface_0_Second.Transposed();
    temp66 = temp6.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
            surface_6_Second((i - 1) * controlPointsRowLength + j, 1) = temp66((i - 1) * controlPointsRowLength + j, 1);
            surface_6_Second((i - 1) * controlPointsRowLength + j, 2) = temp66((i - 1) * controlPointsRowLength + j, 2) + 1;
            surface_6_Second((i - 1) * controlPointsRowLength + j, 3) = temp66((i - 1) * controlPointsRowLength + j, 3) + 1;
        }
    }
    surfaceMatrixVectorSecond.push_back(surface_6_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp6 = F6 * quadMeshMatrix_2.Transposed();
    quadTemp62 = quadTemp6.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp62(i*numberOfQuadMesh+j+1, 2) += 1;
            quadTemp62(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp62);
#endif

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 0 || j == 2) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
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
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 0 || j == 2) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 1 || j == 3) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3);
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 1 || j == 3) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += (eta[j].Value(1)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 1 || j == 3) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3);
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 1 || j == 3) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += (eta[j].Value(2)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 0 || j == 2) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3);
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 0 || j == 2) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += (eta[j].Value(1)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += (eta[j].Value(2)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += eta[j].Value(3);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 1 || j == 3) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3) + 2;
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 1 || j == 3) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += (eta[j].Value(3)+2);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 0 || j == 2) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3) + 2;
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 0 || j == 2) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += (eta[j].Value(1)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += eta[j].Value(2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += (eta[j].Value(3)+2);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 1 || j == 3) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3) + 2;
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 1 || j == 3) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += (eta[j].Value(1)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += (eta[j].Value(2)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += (eta[j].Value(3)+2);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            math_Matrix surfaceTemp(1, (*ctrlpts_1).Size(), 1, 3);
            TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
            math_Matrix tempBC(1, 3, 1, (*ctrlpts_1).Size());
            if (j == 0 || j == 2) {
                tempBC = Q[j] * surfaceMatrixVector[i].Transposed();
            }
            else {
                tempBC = Q[j] * surfaceMatrixVectorSecond[i].Transposed();
            }
            math_Matrix tempBC0 = tempBC.Transposed();
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 1) = tempBC0((p - 1) * controlPointsRowLength + q, 1) + eta[j].Value(1);
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 2) = tempBC0((p - 1) * controlPointsRowLength + q, 2) + eta[j].Value(2) + 2;
                    surfaceTemp((p - 1) * controlPointsRowLength + q, 3) = tempBC0((p - 1) * controlPointsRowLength + q, 3) + eta[j].Value(3) + 2;
                }
            }
            for (Standard_Integer p = 1; p <= controlPointsColLength; p++) {
                for (Standard_Integer q = 1; q <= controlPointsRowLength; q++) {
                    controlPointsTemp.SetValue(p, q, gp_Pnt(surfaceTemp((p - 1) * controlPointsRowLength + q, 1), surfaceTemp((p - 1) * controlPointsRowLength + q, 2), surfaceTemp((p - 1) * controlPointsRowLength + q, 3)));
                }
            }
            controlPointsTotalVector.push_back(scaleControlPoints(controlPointsTemp));

#ifdef USE_QUAD_HEX_MESH
            math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
            if (j == 0 || j == 2) {
                quadTemp_cal1 = Q[j] * quadMeshVector_1[i].Transposed();
            }
            else {
                quadTemp_cal1 = Q[j] * quadMeshVector_2[i].Transposed();
            }
            math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 1) += eta[j].Value(1);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 2) += (eta[j].Value(2)+2);
                    quadTemp_cal2(p*numberOfQuadMesh+q+1, 3) += (eta[j].Value(3)+2);
                }
            }
            for (int p = 0; p < numberOfQuadMesh; p++) {
                for (int q = 0; q < numberOfQuadMesh; q++) {
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 1));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 2));
                    quadMeshPointsTotal.push_back(quadTemp_cal2(p*numberOfQuadMesh+q+1, 3));
                }
            }
#endif

        }
    }

    buildCurvesDiamond2(Uknot, multiplicitiesU, controlPointsColLength, cellSurfacePointsVector, curveControlPointsVector, bSplineCurveVector, controlPointsTotalVector);

    for (int i = 0; i < 384; i++) {
        Handle(Geom_BSplineSurface) bsplineTemp = new Geom_BSplineSurface(
            controlPointsTotalVector[i],
            *Uknot,
            *Vknot,
            (*multiplicitiesU),
            (*multiplicitiesV),
            UDegree,
            VDegree
        );
        cellSurfaceVector.push_back(bsplineTemp);
    }

    GeomFill_BSplineCurves aFillSurface;
    Handle(Geom_BSplineCurve) copyCurveFirst;
    Handle(Geom_BSplineCurve) copyCurveSecond;
    gp_Pnt startPointFirst, startPointSecond, endPointFirst, endPointSecond;
    TColgp_Array1OfPnt curveControl(1, 2);
    TColStd_Array1OfReal curveKnot(1, 2);
    curveKnot.SetValue(1, 0);
    curveKnot.SetValue(2, 1);
    TColStd_Array1OfInteger multiplicitiesCurve(1, 2);
    multiplicitiesCurve.SetValue(1, 2);
    multiplicitiesCurve.SetValue(2, 2);
    Handle(Geom_BSplineCurve) curveFirst, curveSecond;
    std::vector<Handle(Geom_BSplineSurface)> cellBoundarySurfaceVector;
    int edgeBase = 48;

    for (int i = 0; i < 48; i++) {
        copyCurveFirst = Handle(Geom_BSplineCurve)::DownCast(bSplineCurveVector[i]->Copy());
        copyCurveSecond = Handle(Geom_BSplineCurve)::DownCast(bSplineCurveVector[i + edgeBase]->Copy());
        startPointFirst = copyCurveFirst->StartPoint();
        startPointSecond = copyCurveSecond->StartPoint();
        endPointFirst = copyCurveFirst->EndPoint();
        endPointSecond = copyCurveSecond->EndPoint();
        curveControl.SetValue(1, startPointFirst);
        curveControl.SetValue(2, startPointSecond);
        curveFirst = new Geom_BSplineCurve(
            curveControl,
            curveKnot,
            multiplicitiesCurve,
            1
        );
        curveControl.SetValue(1, endPointFirst);
        curveControl.SetValue(2, endPointSecond);
        curveSecond = new Geom_BSplineCurve(
            curveControl,
            curveKnot,
            multiplicitiesCurve,
            1
        );
        aFillSurface = GeomFill_BSplineCurves(copyCurveFirst, copyCurveSecond, curveFirst, curveSecond, GeomFill_CurvedStyle);
        Handle(Geom_BSplineSurface) fillTemp = Handle(Geom_BSplineSurface)::DownCast((aFillSurface.Surface())->Translated(gp_Vec(0, 0, 0)));
        cellBoundarySurfaceVector.push_back(fillTemp);
    }

    for (int i = 0; i < 8; i++) {
        Handle(Geom_BSplineSurface) fillTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Translated(gp_Vec(0, 0.5*cellSizeY, (-0.5)*cellSizeZ)));
        cellBoundarySurfaceVector.push_back(fillTemp);
    }
    for (int i = 8; i < 16; i++) {
        Handle(Geom_BSplineSurface) fillTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Translated(gp_Vec(0, (-0.5)*cellSizeY, (-0.5)*cellSizeZ)));
        cellBoundarySurfaceVector.push_back(fillTemp);
    }
    for (int i = 16; i < 24; i++) {
        Handle(Geom_BSplineSurface) fillTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Translated(gp_Vec((-0.5)*cellSizeX, 0, (-0.5)*cellSizeZ)));
        cellBoundarySurfaceVector.push_back(fillTemp);
    }
    for (int i = 24; i < 32; i++) {
        Handle(Geom_BSplineSurface) fillTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Translated(gp_Vec(0.5*cellSizeX, 0, (-0.5)*cellSizeZ)));
        cellBoundarySurfaceVector.push_back(fillTemp);
    }
    for (int i = 32; i < 40; i++) {
        Handle(Geom_BSplineSurface) fillTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Translated(gp_Vec(0.5*cellSizeX, (-0.5)*cellSizeY, 0)));
        cellBoundarySurfaceVector.push_back(fillTemp);
    }
    for (int i = 40; i < 48; i++) {
        Handle(Geom_BSplineSurface) fillTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Translated(gp_Vec(0.5*cellSizeX, 0.5*cellSizeY, 0)));
        cellBoundarySurfaceVector.push_back(fillTemp);
    }

    std::vector<BRepBuilderAPI_MakeFace> faceMakerVector;
    std::vector<TopoDS_Face> topoFaceVector;
    int cellFaceNumber = 384;
    int internalFaceNumber = cellFaceNumber * this->getCellNumberX() * this->getCellNumberY() * this->getCellNumberZ();
    int coverNumber = 32 * 2 * (this->getCellNumberX() * this->getCellNumberY() + this->getCellNumberX() * this->getCellNumberZ() + this->getCellNumberY() * this->getCellNumberZ());
    int faceNumber = internalFaceNumber + coverNumber;
    std::vector<Handle(Geom_BSplineSurface)> faceVector;

    makeSolidModelFacesDiamond(faceVector, cellSurfaceVector, cellBoundarySurfaceVector);
    
#ifdef USE_QUAD_HEX_MESH
    std::vector<float> quadMeshAll;
    for (int i = 0; i < this->getCellNumberX(); i++) {
        for (int j = 0; j < this->getCellNumberY(); j++) {
            for (int k = 0; k < this->getCellNumberZ(); k++) {
                for (int m = 0; m < numberOfQuadMesh*numberOfQuadMesh*384; m++) {
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3]+i*cellSizeX);
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3+1]+j*cellSizeY);
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3+2]+k*cellSizeZ);
                }
            }
        }
    }
    writeHexahedralMesh(quadMeshAll, 384*this->getCellNumberX()*this->getCellNumberY()*this->getCellNumberZ(), 192, numberOfQuadMesh, "diamond_hex.off");

    std::vector<float> quadMeshEdge;
    for (int i=0;i<faceVector.size();i++) {
        quadMeshEdge.clear();
        extractMesh(faceVector[i], quadMeshEdge, numberOfQuadMesh);
        for (int j = 0; j < numberOfQuadMesh*numberOfQuadMesh*3; j++) {
            quadMeshAll.push_back(quadMeshEdge[j]);
        }
    }
    writeMesh(quadMeshAll, 384*this->getCellNumberX()*this->getCellNumberY()*this->getCellNumberZ()+faceVector.size(), numberOfQuadMesh, "diamond.off");
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
                faceVector[i]->SetPole(j,k,gp_Pnt((faceVector[i]->Pole(j,k)).X()*(1+scaleFactor*(faceVector[i]->Pole(j,k)).X()),
                (faceVector[i]->Pole(j,k)).Y(), (faceVector[i]->Pole(j,k)).Z()));
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
        std::cerr << "TPMS2STEP > Failed to create Diamond TPMS solid model." << std::endl;
        return false;
    }
    cout << "TPMS2STEP > Diamond TPMS solid model generated successfully." << endl;
    return true;
}

void TPMS_D::errorControlledSampling(vector<float>& knotTemp, float* knot, int* NumRows, int* NumCols, int* numberOfPointCloud, Eigen::MatrixXf& surfacePlusMatrix,
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
        cout << "TPMS2STEP > Diamond Error Control Iteration: " << iterationCount << endl;
        int length = (paramRuler.size()+1)/2;// the number of control points is the same as the number of sample points
        int halfLength = (length+1)/2;
        sampleDiamond(paramRuler);
        pythonEXEC(this->modelType, this->getOffsetValue1(), this->getOffsetValue2());
        (*NumRows) = (paramRuler.size()+1)/2;
        (*NumCols) = (paramRuler.size()+1)/2;
        controlPointsPlusPIA.resize((*NumRows)*(*NumCols), 3);
        controlPointsMinusPIA.resize((*NumRows)*(*NumCols), 3);
        (*numberOfPointCloud) = (paramRuler.size()+1)/2*(paramRuler.size()+1)/2;
        U.resize(*numberOfPointCloud);
        V.resize(*numberOfPointCloud);
        surfacePlusMatrix.resize((*numberOfPointCloud), 3);
        surfaceMinusMatrix.resize((*numberOfPointCloud), 3);
        // read sample points
        readSemiFundamentalPatch(surfacePlusMatrix, surfaceMinusMatrix, numberOfPointCloud, this->getOffsetValue1());
        ofs.open("../data/uv/uv_Diamond_pia.off", ios::out);
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
            knot = (float*)malloc(((paramRuler.size()+1)/2+4)*sizeof(float));
            knotTemp.clear();
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

void TPMS_D::readSemiFundamentalPatch(Eigen::MatrixXf& surfacePlusMatrix, Eigen::MatrixXf& surfaceMinusMatrix, int* numberOfPointCloud, float offsetValue) {
    string tempString;
    ifstream ifs;
    int vertex, edge, face;
    float a, b, c;
    ifs = ifstream("../data/SFP/Diamond/SFP_Diamond_offset_" + to_string(int(offsetValue * (-100) + 1 - 1)) + ".off", ios::in);
    ifs >> tempString;
    ifs >> vertex >> edge >> face;
    (*numberOfPointCloud) = vertex;
    for (int i = 0; i < vertex; i++) {
        ifs >> a >> b >> c;
        surfacePlusMatrix(i, 0) = a;
        surfacePlusMatrix(i, 1) = b;
        surfacePlusMatrix(i, 2) = c;
    }
    ifs.close();
    ifs = ifstream("../data/SFP/Diamond/SFP_Diamond_offset_" + to_string(int((offsetValue * (-100) + 1 - 1) * (-1))) + ".off", ios::in);
    ifs >> tempString;
    ifs >> vertex >> edge >> face;
    for (int i = 0; i < vertex; i++) {
        ifs >> a >> b >> c;
        surfaceMinusMatrix(i, 0) = a;
        surfaceMinusMatrix(i, 1) = b;
        surfaceMinusMatrix(i, 2) = c;
    }
    ifs.close();
}


void TPMS_D::buildCurvesDiamond(TColStd_Array1OfReal* Uknot, TColStd_Array1OfInteger* multiplicitiesU, int controlPointsColLength, std::vector<TColgp_Array2OfPnt>& cellSurfacePointsVector,
std::vector<TColgp_Array1OfPnt>& curveControlPointsVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector, std::vector<TColgp_Array2OfPnt>& controlPointsTotalVector) {
    // x axis
    TColgp_Array1OfPnt curveCpts1(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[107].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts1);
    Handle(Geom_BSplineCurve) curve1 = new Geom_BSplineCurve(
        curveCpts1,
        (*Uknot),
        (*multiplicitiesU),
        3
    );
    TColgp_Array1OfPnt curveCpts2(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[119].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    Handle(Geom_BSplineCurve) curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts3(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[119].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    Handle(Geom_BSplineCurve) curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts4(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[111].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts4);
    Handle(Geom_BSplineCurve) curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts5(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[98].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    Handle(Geom_BSplineCurve) curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts6(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[102].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    Handle(Geom_BSplineCurve) curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts7(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[102].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    Handle(Geom_BSplineCurve) curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts8(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[114].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts8);
    Handle(Geom_BSplineCurve) curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );
    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[170].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[174].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[174].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[186].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[179].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[191].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[191].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[183].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    // y axis
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[122].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[130].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[130].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[142].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[125].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[137].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[137].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[133].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[98].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[106].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[106].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[118].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[101].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[113].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[113].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[109].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    // z axis
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[59].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[51].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[51].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[55].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[53].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[49].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[49].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[57].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[9].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[1].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[1].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[5].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[7].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[3].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[3].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[11].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );
    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);
}

void TPMS_D::buildCurvesDiamond2(TColStd_Array1OfReal* Uknot, TColStd_Array1OfInteger* multiplicitiesU, int controlPointsColLength, std::vector<TColgp_Array2OfPnt>& cellSurfacePointsVector,
std::vector<TColgp_Array1OfPnt>& curveControlPointsVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector, std::vector<TColgp_Array2OfPnt>& controlPointsTotalVector) {
    TColgp_Array1OfPnt curveCpts1(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts2(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts3(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts4(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts5(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts6(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts7(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts8(1, controlPointsColLength);

    Handle(Geom_BSplineCurve) curve1, curve2, curve3, curve4, curve5, curve6, curve7, curve8, curve9, curve10, curve11, curve12, curve13, curve14, curve15, curve16, curve17, curve18;
    // x axis
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[299].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[311].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[311].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[303].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[290].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[294].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[294].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[306].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[362].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[366].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[366].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[378].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[371].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[383].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[383].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[375].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    // y axis
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[314].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[322].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[322].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[334].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[317].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[329].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[329].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[325].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[290].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[298].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[298].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[310].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[293].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[305].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[305].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[301].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    // z axis
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[251].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[243].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[243].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[247].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[245].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[241].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[241].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[249].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, controlPointsTotalVector[201].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts1);
    curve1 = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, controlPointsTotalVector[193].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    curve2 = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, controlPointsTotalVector[193].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts3);
    curve3 = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, controlPointsTotalVector[197].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    curve4 = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, controlPointsTotalVector[199].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    curve5 = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, controlPointsTotalVector[195].Value(controlPointsColLength, i));
    }
    curveControlPointsVector.push_back(curveCpts6);
    curve6 = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, controlPointsTotalVector[195].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts7);
    curve7 = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, controlPointsTotalVector[203].Value(i, controlPointsColLength));
    }
    curveControlPointsVector.push_back(curveCpts8);
    curve8 = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curve1);
    bSplineCurveVector.push_back(curve2);
    bSplineCurveVector.push_back(curve3);
    bSplineCurveVector.push_back(curve4);
    bSplineCurveVector.push_back(curve5);
    bSplineCurveVector.push_back(curve6);
    bSplineCurveVector.push_back(curve7);
    bSplineCurveVector.push_back(curve8);
}

void TPMS_D::makeSolidModelFacesDiamond(std::vector<Handle(Geom_BSplineSurface)>& faceVector, std::vector<Handle(Geom_BSplineSurface)>& cellSurfaceVector, std::vector<Handle(Geom_BSplineSurface)>& cellBoundarySurfaceVector) {
    for (int i = 0; i < 16; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Copy());
        for (int i = 0; i < this->getCellNumberY(); i++) {
            for (int j = 0; j < this->getCellNumberZ(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(0, i * cellSizeY, j * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(this->getCellNumberX() * cellSizeX, 0, 0)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }
    for (int i = 48; i < 64; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Copy());
        for (int i = 0; i < this->getCellNumberY(); i++) {
            for (int j = 0; j < this->getCellNumberZ(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(0, i * cellSizeY, j * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(this->getCellNumberX() * cellSizeX, 0, 0)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }

    for (int i = 16; i < 32; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Copy());
        for (int i = 0; i < this->getCellNumberX(); i++) {
            for (int j = 0; j < this->getCellNumberZ(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(i * cellSizeX, 0, j * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, this->getCellNumberY() * cellSizeY, 0)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }

    for (int i = 64; i < 80; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Copy());
        for (int i = 0; i < this->getCellNumberX(); i++) {
            for (int j = 0; j < this->getCellNumberZ(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(i * cellSizeX, 0, j * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, this->getCellNumberY() * cellSizeY, 0)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }

    for (int i = 32; i < 48; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Copy());
        for (int i = 0; i < this->getCellNumberX(); i++) {
            for (int j = 0; j < this->getCellNumberY(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(i * cellSizeX, j * cellSizeY, 0)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, 0, this->getCellNumberZ() * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }

    for (int i = 80; i < 96; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i]->Copy());
        for (int i = 0; i < this->getCellNumberX(); i++) {
            for (int j = 0; j < this->getCellNumberY(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(i * cellSizeX, j * cellSizeY, 0)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, 0, this->getCellNumberZ() * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }
}