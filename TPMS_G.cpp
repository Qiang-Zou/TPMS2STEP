#include "TPMS_G.h"

TPMS_G::TPMS_G()
{
}

TPMS_G::~TPMS_G()
{
}

bool TPMS_G::constructSolidModel() {
    // define some params
    this->maxIteration = 1;
    this->setOffsetValue(0.1, 0);
    this->setCellNumbers(1, 1, 1);
    this->setTolerance(0.005);
    this->cellScaleX = 1.0;
    this->cellScaleY = 1.0;
    this->cellScaleZ = 1.0;
    float epsilon_1 = 0.5*this->getTolerance();
    float epsilon_2 = 0.5*this->getTolerance();

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
    this->cellSizeX = 2.0 * this->cellScaleX;
    this->cellSizeY = 2.0 * this->cellScaleY;
    this->cellSizeZ = 2.0 * this->cellScaleZ;
    this->minControlPointNumberPerEdge = 10;
    float* knot;

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
    // error-controlled sampling
    errorControlledSampling(knotTemp, knot, &NumRows, &NumCols, &numberOfPointCloud, surfacePlusMatrix, surfaceMinusMatrix, U, V, controlPointsPlusPIA, controlPointsMinusPIA, n, epsilon_2, paramRuler);
    cout << "TPMS2STEP > Error Control ok." << endl;
    if (knot) {
        free(knot);
    }
    knot = (float*)malloc(knotTemp.size()*sizeof(float));
    for (int i=0;i<knotTemp.size();i++) {
        knot[i] = knotTemp[i];
    }
    // assemble
    // define knot vectors
    TColStd_Array1OfReal* Uknot = new TColStd_Array1OfReal(1, paramRuler.size()-2);
    TColStd_Array1OfReal* Vknot = new TColStd_Array1OfReal(1, paramRuler.size()-2);

    for (int i=1;i<=paramRuler.size()-2;i++) {
        (*Uknot).SetValue(i,(i-1)*1.0/(paramRuler.size()-3));
        (*Vknot).SetValue(i,(i-1)*1.0/(paramRuler.size()-3));
    }
    TColStd_Array1OfInteger* multiplicitiesU = new TColStd_Array1OfInteger(1, paramRuler.size()-2);
    TColStd_Array1OfInteger* multiplicitiesV = new TColStd_Array1OfInteger(1, paramRuler.size()-2);
    for (int i=2;i<=paramRuler.size()-3;i++) {
        (*multiplicitiesU).SetValue(i,1);
        (*multiplicitiesV).SetValue(i,1);
    }
    (*multiplicitiesU).SetValue(1, 4);
    (*multiplicitiesU).SetValue(paramRuler.size()-2, 4);
    (*multiplicitiesV).SetValue(1, 4);
    (*multiplicitiesV).SetValue(paramRuler.size()-2, 4);

#ifdef USE_QUAD_HEX_MESH
    int numberOfQuadMesh = 12;
    std::vector<float> quadMeshPoints_1;
    std::vector<float> quadMeshPoints_2;
    std::vector<float> quadMeshPointsTotal;
#endif

    std::vector<Handle(Geom_BSplineSurface)> cellSurfaceVector;
    std::vector<TColgp_Array2OfPnt> controlPointsVector;
    std::vector<TColgp_Array2OfPnt> cellSurfacePointsVector;

    TColgp_Array2OfPnt* ctrlpts_1 = new TColgp_Array2OfPnt(1, NumRows, 1, NumCols);
    TColgp_Array2OfPnt* ctrlpts_2 = new TColgp_Array2OfPnt(1, NumRows, 1, NumCols);
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
    controlPointsVector.push_back(scaleControlPoints((*ctrlpts_1)));
    cellSurfacePointsVector.push_back(scaleControlPoints((*ctrlpts_1)));

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

    int controlPointsRowLength = NumRows;
    int controlPointsColLength = NumCols;

    std::vector<math_Matrix> surfaceMatrixVector;
    std::vector<math_Matrix> surfaceMatrixSecondVector;
    math_Matrix surface_1(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_2(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_3(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_4(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_5(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_6(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_1_Second(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_2_Second(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_3_Second(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_4_Second(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_5_Second(1, NumRows*NumCols, 1, 3);
    math_Matrix surface_6_Second(1, NumRows*NumCols, 1, 3);

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
    quadMeshVector_1.push_back(quadMeshMatrix_1);
    quadMeshVector_2.push_back(quadMeshMatrix_2);

    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,1));
            quadMeshPointsTotal.push_back(quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,2));
            quadMeshPointsTotal.push_back(quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,3));
        }
    }
#endif

    TColgp_Array2OfPnt controlPoints_2(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_3(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_4(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_5(1, controlPointsColLength, 1, controlPointsRowLength);
    TColgp_Array2OfPnt controlPoints_6(1, controlPointsColLength, 1, controlPointsRowLength);

    math_Matrix Tx(1, 3, 1, 3, 0);
    math_Matrix Ty(1, 3, 1, 3, 0);
    math_Matrix Tz(1, 3, 1, 3, 0);
    Tx(1, 1) = -1; Tx(2, 2) = 1; Tx(3, 3) = 1;
    Ty(1, 1) = 1; Ty(2, 2) = -1; Ty(3, 3) = 1;
    Tz(1, 1) = 1; Tz(2, 2) = 1; Tz(3, 3) = -1;
    math_Matrix Txy = Tx * Ty;
    math_Matrix Txz = Tx * Tz;
    math_Matrix Tyz = Ty * Tz;

    math_Matrix D = Tz * Ry((-0.5) * M_PI) * Rx((-0.5) * M_PI);

    for (int i = 1; i <= NumCols; i++) {
        for (int j = 1; j <= NumRows; j++) {
            surface_1((i - 1) * NumRows + j, 1) = (*ctrlpts_1).Value(i, j).X();
            surface_1((i - 1) * NumRows + j, 2) = (*ctrlpts_1).Value(i, j).Y();
            surface_1((i - 1) * NumRows + j, 3) = (*ctrlpts_1).Value(i, j).Z();
        }
    }
    surfaceMatrixVector.push_back(surface_1);

    for (int i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (int j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_1_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1) = (*ctrlpts_2).Value(i, j).X();
            surface_1_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2) = (*ctrlpts_2).Value(i, j).Y();
            surface_1_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) = (*ctrlpts_2).Value(i, j).Z();
        }
    }
    surfaceMatrixSecondVector.push_back(surface_1_Second);

    math_Matrix temp12 = D * surface_1.Transposed();
    math_Matrix temp122 = temp12.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_2((i - 1) * controlPointsColLength + j, 1) = temp122((i - 1) * controlPointsColLength + j, 1);
            surface_2((i - 1) * controlPointsColLength + j, 2) = temp122((i - 1) * controlPointsColLength + j, 2);
            surface_2((i - 1) * controlPointsColLength + j, 3) = temp122((i - 1) * controlPointsColLength + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_2);

    math_Matrix temp12_Second = D * surface_1_Second.Transposed();
    math_Matrix temp122_Second = temp12_Second.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsColLength; j++) {
            surface_2_Second((i - 1) * controlPointsColLength + j, 1) = temp122_Second((i - 1) * controlPointsColLength + j, 1);
            surface_2_Second((i - 1) * controlPointsColLength + j, 2) = temp122_Second((i - 1) * controlPointsColLength + j, 2);
            surface_2_Second((i - 1) * controlPointsColLength + j, 3) = temp122_Second((i - 1) * controlPointsColLength + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_2_Second);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp2 = D * quadMeshMatrix_1.Transposed();
    math_Matrix quadTemp22 = quadTemp2.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp22(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp22);
    quadTemp2 = D * quadMeshMatrix_2.Transposed();
    quadTemp22 = quadTemp2.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp22(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp22);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp22(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp22(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp22(i*numberOfQuadMesh+j+1, 3));
        }
    }
#endif

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsColLength; j++) {
            controlPoints_2.SetValue(i, j, gp_Pnt(surface_2_Second((i - 1) * controlPointsColLength + j, 1), surface_2_Second((i - 1) * controlPointsColLength + j, 2), surface_2_Second((i - 1) * controlPointsColLength + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_2));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_2));

    math_Matrix temp23 = D * surface_2.Transposed();
    math_Matrix temp233 = temp23.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_3((i - 1) * NumRows + j, 1) = temp233((i - 1) * NumRows + j, 1);
            surface_3((i - 1) * NumRows + j, 2) = temp233((i - 1) * NumRows + j, 2);
            surface_3((i - 1) * NumRows + j, 3) = temp233((i - 1) * NumRows + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_3);

    math_Matrix temp23_Second = D * surface_2_Second.Transposed();
    math_Matrix temp233_Second = temp23_Second.Transposed();
    for (Standard_Integer i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_3_Second((i - 1) * controlPointsColLength + j, 1) = temp233_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1);
            surface_3_Second((i - 1) * controlPointsColLength + j, 2) = temp233_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2);
            surface_3_Second((i - 1) * controlPointsColLength + j, 3) = temp233_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_3_Second);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp3 = D * quadMeshVector_1[quadMeshVector_1.size()-1].Transposed();
    math_Matrix quadTemp32 = quadTemp3.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp32(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp32);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp32(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp32(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp32(i*numberOfQuadMesh+j+1, 3));
        }
    }
    quadTemp3 = D * quadMeshVector_2[quadMeshVector_2.size()-1].Transposed();
    quadTemp32 = quadTemp3.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp32(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp32);
#endif

    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            controlPoints_3.SetValue(i, j, gp_Pnt(surface_3((i - 1) * NumRows + j, 1), surface_3((i - 1) * NumRows + j, 2), surface_3((i - 1) * NumRows + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_3));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_3));

    math_Matrix temp34 = D * surface_3.Transposed();
    math_Matrix temp344 = temp34.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_4((i - 1) * NumRows + j, 1) = temp344((i - 1) * NumRows + j, 1);
            surface_4((i - 1) * NumRows + j, 2) = temp344((i - 1) * NumRows + j, 2);
            surface_4((i - 1) * NumRows + j, 3) = temp344((i - 1) * NumRows + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_4);

    math_Matrix temp34_Second = D * surface_3_Second.Transposed();
    math_Matrix temp344_Second = temp34_Second.Transposed();
    for (Standard_Integer i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_4_Second((i - 1) * controlPointsColLength + j, 1) = temp344_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1);
            surface_4_Second((i - 1) * controlPointsColLength + j, 2) = temp344_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2);
            surface_4_Second((i - 1) * controlPointsColLength + j, 3) = temp344_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_4_Second);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp4 = D * quadMeshVector_1[quadMeshVector_1.size()-1].Transposed();
    math_Matrix quadTemp42 = quadTemp4.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp42(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp42);
    quadTemp4 = D * quadMeshVector_2[quadMeshVector_2.size()-1].Transposed();
    quadTemp42 = quadTemp4.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp42(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp42);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp42(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp42(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp42(i*numberOfQuadMesh+j+1, 3));
        }
    }
#endif

    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            controlPoints_4.SetValue(i, j, gp_Pnt(surface_4_Second((i - 1) * NumRows + j, 1), surface_4_Second((i - 1) * NumRows + j, 2), surface_4_Second((i - 1) * NumRows + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_4));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_4));

    math_Matrix temp45 = D * surface_4.Transposed();
    math_Matrix temp455 = temp45.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_5((i - 1) * NumRows + j, 1) = temp455((i - 1) * NumRows + j, 1);
            surface_5((i - 1) * NumRows + j, 2) = temp455((i - 1) * NumRows + j, 2);
            surface_5((i - 1) * NumRows + j, 3) = temp455((i - 1) * NumRows + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_5);

    math_Matrix temp45_Second = D * surface_4_Second.Transposed();
    math_Matrix temp455_Second = temp45_Second.Transposed();
    for (Standard_Integer i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_5_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1) = temp455_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1);
            surface_5_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2) = temp455_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2);
            surface_5_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) = temp455_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_5_Second);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp5 = D * quadMeshVector_1[quadMeshVector_1.size()-1].Transposed();
    math_Matrix quadTemp52 = quadTemp5.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp52(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp52);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp52(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp52(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp52(i*numberOfQuadMesh+j+1, 3));
        }
    }
    quadTemp5 = D * quadMeshVector_2[quadMeshVector_2.size()-1].Transposed();
    quadTemp52 = quadTemp5.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp52(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp52);
#endif

    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            controlPoints_5.SetValue(i, j, gp_Pnt(surface_5((i - 1) * NumRows + j, 1), surface_5((i - 1) * NumRows + j, 2), surface_5((i - 1) * NumRows + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_5));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_5));

    math_Matrix temp56 = D * surface_5.Transposed();
    math_Matrix temp566 = temp56.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_6((i - 1) * NumRows + j, 1) = temp566((i - 1) * NumRows + j, 1);
            surface_6((i - 1) * NumRows + j, 2) = temp566((i - 1) * NumRows + j, 2);
            surface_6((i - 1) * NumRows + j, 3) = temp566((i - 1) * NumRows + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_6);

    math_Matrix temp56_Second = D * surface_5_Second.Transposed();
    math_Matrix temp566_Second = temp56_Second.Transposed();
    for (Standard_Integer i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_6_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1) = temp566_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1);
            surface_6_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2) = temp566_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2);
            surface_6_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) = temp566_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_6_Second);

#ifdef USE_QUAD_HEX_MESH
    math_Matrix quadTemp6 = D * quadMeshVector_1[quadMeshVector_1.size()-1].Transposed();
    math_Matrix quadTemp62 = quadTemp6.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp62(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp62);
    quadTemp6 = D * quadMeshVector_2[quadMeshVector_2.size()-1].Transposed();
    quadTemp62 = quadTemp6.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp62(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp62);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp62(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp62(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp62(i*numberOfQuadMesh+j+1, 3));
        }
    }
#endif

    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            controlPoints_6.SetValue(i, j, gp_Pnt(surface_6_Second((i - 1) * NumRows + j, 1), surface_6_Second((i - 1) * NumRows + j, 2), surface_6_Second((i - 1) * NumRows + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_6));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_6));

    // copy the bounding cell 3 times
    for (int i = 0; i < 6; i++) {
        math_Matrix surfaceTemp(1, 3, 1, NumRows*NumCols);
        if (i % 2 == 0) {
            surfaceTemp = Ty * surfaceMatrixSecondVector[i].Transposed();
        }
        else {
            surfaceTemp = Ty * surfaceMatrixVector[i].Transposed();
        }
        math_Matrix surfaceTemp2 = surfaceTemp.Transposed();
        math_Matrix surfaceTemp3(1, controlPointsColLength*controlPointsRowLength, 1, 3);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 1) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 1) + 1;
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 2) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 2);
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 3) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 3);
            }
        }
        surfaceMatrixVector.push_back(surfaceTemp3);

#ifdef USE_QUAD_HEX_MESH
        math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
        if (i % 2 == 0) {
            quadTemp_cal1 = Ty * quadMeshVector_2[i].Transposed();
        }
        else {
            quadTemp_cal1 = Ty * quadMeshVector_1[i].Transposed();
        }
        math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 1) += 1;
            }
        }
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 1));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 2));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 3));
            }
        }
#endif

        TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                controlPointsTemp.SetValue(i, j, gp_Pnt(surfaceTemp3((i - 1) * controlPointsRowLength + j, 1), surfaceTemp3((i - 1) * controlPointsRowLength + j, 2), surfaceTemp3((i - 1) * controlPointsRowLength + j, 3)));
            }
        }
        controlPointsVector.push_back(scaleControlPoints(controlPointsTemp));
        cellSurfacePointsVector.push_back(scaleControlPoints(controlPointsTemp));
    }
    

    for (int i = 0; i < 6; i++) {
        math_Matrix surfaceTemp(1, 3, 1, NumRows*NumCols);
        if (i % 2 == 1) {
            surfaceTemp = Tyz * surfaceMatrixSecondVector[i].Transposed();
        }
        else {
            surfaceTemp = Tyz * surfaceMatrixVector[i].Transposed();
        }
        math_Matrix surfaceTemp2 = surfaceTemp.Transposed();
        math_Matrix surfaceTemp3(1, controlPointsColLength * controlPointsRowLength, 1, 3);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 1) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 1) + 1;
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 2) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 2) + 1;
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 3) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 3);
            }
        }
        surfaceMatrixVector.push_back(surfaceTemp3);

#ifdef USE_QUAD_HEX_MESH
        math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
        if (i % 2 == 1) {
            quadTemp_cal1 = Tyz * quadMeshVector_2[i].Transposed();
        }
        else {
            quadTemp_cal1 = Tyz * quadMeshVector_1[i].Transposed();
        }
        math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 1) += 1;
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 2) += 1;
            }
        }
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 1));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 2));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 3));
            }
        }
#endif

        TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                controlPointsTemp.SetValue(i, j, gp_Pnt(surfaceTemp3((i - 1) * controlPointsRowLength + j, 1), surfaceTemp3((i - 1) * controlPointsRowLength + j, 2), surfaceTemp3((i - 1) * controlPointsRowLength + j, 3)));
            }
        }
        controlPointsVector.push_back(scaleControlPoints(controlPointsTemp));
        cellSurfacePointsVector.push_back(scaleControlPoints(controlPointsTemp));
    }

    for (int i = 0; i < 6; i++) {
        math_Matrix surfaceTemp(1, 3, 1, NumRows*NumCols);
        if (i % 2 == 1) {
            surfaceTemp = Txy * surfaceMatrixSecondVector[i].Transposed();
        }
        else {
            surfaceTemp = Txy * surfaceMatrixVector[i].Transposed();
        }
        math_Matrix surfaceTemp2 = surfaceTemp.Transposed();
        math_Matrix surfaceTemp3(1, controlPointsColLength * controlPointsRowLength, 1, 3);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 1) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 1) + 1;
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 2) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 2);
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 3) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 3) - 1;
            }
        }
        surfaceMatrixVector.push_back(surfaceTemp3);

#ifdef USE_QUAD_HEX_MESH
        math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
        if (i % 2 == 1) {
            quadTemp_cal1 = Txy * quadMeshVector_2[i].Transposed();
        }
        else {
            quadTemp_cal1 = Txy * quadMeshVector_1[i].Transposed();
        }
        math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 1) += 1;
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 3) -= 1;
            }
        }
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 1));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 2));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 3));
            }
        }
#endif

        TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                controlPointsTemp.SetValue(i, j, gp_Pnt(surfaceTemp3((i - 1) * controlPointsRowLength + j, 1), surfaceTemp3((i - 1) * controlPointsRowLength + j, 2), surfaceTemp3((i - 1) * controlPointsRowLength + j, 3)));
            }
        }
        controlPointsVector.push_back(scaleControlPoints(controlPointsTemp));
        cellSurfacePointsVector.push_back(scaleControlPoints(controlPointsTemp));
    }

    std::vector<TColgp_Array1OfPnt> curveControlPointsVector;
    std::vector<Handle(Geom_BSplineCurve)> bSplineCurveVector;
    buildCurvesGyroid(Uknot, multiplicitiesU, controlPointsColLength, cellSurfacePointsVector, curveControlPointsVector, bSplineCurveVector);
    surfaceMatrixVector.clear();
    surfaceMatrixSecondVector.clear();
    // read the sample points
    readSemiFundamentalPatch(surfacePlusMatrix, surfaceMinusMatrix, &numberOfPointCloud, this->getOffsetValue2());
    constrainedPIAGPU(knot, 3, numberOfPointCloud, NumRows, NumCols, surfacePlusMatrix, surfaceMinusMatrix, U, V,
    &controlPointsPlusPIA, &controlPointsMinusPIA, this->getOffsetValue2(), this->modelType);
    for (int i=0;i<NumRows*NumCols;i++) {
        controlPointsPlusPIA(i,0) += 1;
        controlPointsPlusPIA(i,1) += 0.5;
        controlPointsMinusPIA(i,0) += 1;
        controlPointsMinusPIA(i,1) += 0.5;
    }
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
    if (knot) {
        free(knot);
    }
    controlPointsVector.push_back(scaleControlPoints((*ctrlpts_1)));
    cellSurfacePointsVector.push_back(scaleControlPoints((*ctrlpts_1)));

#ifdef USE_QUAD_HEX_MESH
    quadMeshPoints_1.clear();
    quadMeshPoints_2.clear();
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
    quadMeshVector_1.clear();
    quadMeshVector_2.clear();

    quadMeshVector_1.push_back(quadMeshMatrix_1);
    quadMeshVector_2.push_back(quadMeshMatrix_2);

    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,1));
            quadMeshPointsTotal.push_back(quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,2));
            quadMeshPointsTotal.push_back(quadMeshMatrix_1(i*numberOfQuadMesh+j+1 ,3));
        }
    }
#endif
    controlPointsColLength = NumCols;
    controlPointsRowLength = NumRows;

    surface_1 = math_Matrix(1, NumRows*NumCols, 1, 3);
    surface_2 = math_Matrix(1, NumRows*NumCols, 1, 3);
    surface_3 = math_Matrix(1, NumRows*NumCols, 1, 3);
    surface_4 = math_Matrix(1, NumRows*NumCols, 1, 3);
    surface_5 = math_Matrix(1, NumRows*NumCols, 1, 3);
    surface_6 = math_Matrix(1, NumRows*NumCols, 1, 3);

    controlPoints_2 = TColgp_Array2OfPnt(1, controlPointsColLength, 1, controlPointsRowLength);
    controlPoints_3 = TColgp_Array2OfPnt(1, controlPointsColLength, 1, controlPointsRowLength);
    controlPoints_4 = TColgp_Array2OfPnt(1, controlPointsColLength, 1, controlPointsRowLength);
    controlPoints_5 = TColgp_Array2OfPnt(1, controlPointsColLength, 1, controlPointsRowLength);
    controlPoints_6 = TColgp_Array2OfPnt(1, controlPointsColLength, 1, controlPointsRowLength);
    for (int i = 1; i <= NumCols; i++) {
        for (int j = 1; j <= NumRows; j++) {
            surface_1((i - 1) * NumRows + j, 1) = (*ctrlpts_1).Value(i, j).X();
            surface_1((i - 1) * NumRows + j, 2) = (*ctrlpts_1).Value(i, j).Y();
            surface_1((i - 1) * NumRows + j, 3) = (*ctrlpts_1).Value(i, j).Z();
        }
    }
    surfaceMatrixVector.push_back(surface_1);

    for (int i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (int j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_1_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1) = (*ctrlpts_2).Value(i, j).X();
            surface_1_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2) = (*ctrlpts_2).Value(i, j).Y();
            surface_1_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) = (*ctrlpts_2).Value(i, j).Z();
        }
    }
    surfaceMatrixSecondVector.push_back(surface_1_Second);

    temp12 = D * surface_1.Transposed();
    temp122 = temp12.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_2((i - 1) * controlPointsColLength + j, 1) = temp122((i - 1) * controlPointsColLength + j, 1);
            surface_2((i - 1) * controlPointsColLength + j, 2) = temp122((i - 1) * controlPointsColLength + j, 2);
            surface_2((i - 1) * controlPointsColLength + j, 3) = temp122((i - 1) * controlPointsColLength + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_2);

    temp12_Second = D * surface_1_Second.Transposed();
    temp122_Second = temp12_Second.Transposed();
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsColLength; j++) {
            surface_2_Second((i - 1) * controlPointsColLength + j, 1) = temp122_Second((i - 1) * controlPointsColLength + j, 1);
            surface_2_Second((i - 1) * controlPointsColLength + j, 2) = temp122_Second((i - 1) * controlPointsColLength + j, 2);
            surface_2_Second((i - 1) * controlPointsColLength + j, 3) = temp122_Second((i - 1) * controlPointsColLength + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_2_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp2 = D * quadMeshMatrix_1.Transposed();
    quadTemp22 = quadTemp2.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp22(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp22);
    quadTemp2 = D * quadMeshMatrix_2.Transposed();
    quadTemp22 = quadTemp2.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp22(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp22);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp22(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp22(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp22(i*numberOfQuadMesh+j+1, 3));
        }
    }
#endif

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        for (Standard_Integer j = 1; j <= controlPointsColLength; j++) {
            controlPoints_2.SetValue(i, j, gp_Pnt(surface_2_Second((i - 1) * controlPointsColLength + j, 1), surface_2_Second((i - 1) * controlPointsColLength + j, 2), surface_2_Second((i - 1) * controlPointsColLength + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_2));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_2));

    temp23 = D * surface_2.Transposed();
    temp233 = temp23.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_3((i - 1) * NumRows + j, 1) = temp233((i - 1) * NumRows + j, 1);
            surface_3((i - 1) * NumRows + j, 2) = temp233((i - 1) * NumRows + j, 2);
            surface_3((i - 1) * NumRows + j, 3) = temp233((i - 1) * NumRows + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_3);

    temp23_Second = D * surface_2_Second.Transposed();
    temp233_Second = temp23_Second.Transposed();
    for (Standard_Integer i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_3_Second((i - 1) * controlPointsColLength + j, 1) = temp233_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1);
            surface_3_Second((i - 1) * controlPointsColLength + j, 2) = temp233_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2);
            surface_3_Second((i - 1) * controlPointsColLength + j, 3) = temp233_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_3_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp3 = D * quadMeshVector_1[quadMeshVector_1.size()-1].Transposed();
    quadTemp32 = quadTemp3.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp32(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp32);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp32(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp32(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp32(i*numberOfQuadMesh+j+1, 3));
        }
    }
    quadTemp3 = D * quadMeshVector_2[quadMeshVector_2.size()-1].Transposed();
    quadTemp32 = quadTemp3.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp32(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp32);
#endif
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            controlPoints_3.SetValue(i, j, gp_Pnt(surface_3((i - 1) * NumRows + j, 1), surface_3((i - 1) * NumRows + j, 2), surface_3((i - 1) * NumRows + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_3));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_3));

    temp34 = D * surface_3.Transposed();
    temp344 = temp34.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_4((i - 1) * NumRows + j, 1) = temp344((i - 1) * NumRows + j, 1);
            surface_4((i - 1) * NumRows + j, 2) = temp344((i - 1) * NumRows + j, 2);
            surface_4((i - 1) * NumRows + j, 3) = temp344((i - 1) * NumRows + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_4);

    temp34_Second = D * surface_3_Second.Transposed();
    temp344_Second = temp34_Second.Transposed();
    for (Standard_Integer i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_4_Second((i - 1) * controlPointsColLength + j, 1) = temp344_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1);
            surface_4_Second((i - 1) * controlPointsColLength + j, 2) = temp344_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2);
            surface_4_Second((i - 1) * controlPointsColLength + j, 3) = temp344_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_4_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp4 = D * quadMeshVector_1[quadMeshVector_1.size()-1].Transposed();
    quadTemp42 = quadTemp4.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp42(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp42);
    quadTemp4 = D * quadMeshVector_2[quadMeshVector_2.size()-1].Transposed();
    quadTemp42 = quadTemp4.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp42(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp42);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp42(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp42(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp42(i*numberOfQuadMesh+j+1, 3));
        }
    }
#endif
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            controlPoints_4.SetValue(i, j, gp_Pnt(surface_4_Second((i - 1) * NumRows + j, 1), surface_4_Second((i - 1) * NumRows + j, 2), surface_4_Second((i - 1) * NumRows + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_4));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_4));

    temp45 = D * surface_4.Transposed();
    temp455 = temp45.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_5((i - 1) * NumRows + j, 1) = temp455((i - 1) * NumRows + j, 1);
            surface_5((i - 1) * NumRows + j, 2) = temp455((i - 1) * NumRows + j, 2);
            surface_5((i - 1) * NumRows + j, 3) = temp455((i - 1) * NumRows + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_5);

    temp45_Second = D * surface_4_Second.Transposed();
    temp455_Second = temp45_Second.Transposed();
    for (Standard_Integer i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_5_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1) = temp455_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1);
            surface_5_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2) = temp455_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2);
            surface_5_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) = temp455_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_5_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp5 = D * quadMeshVector_1[quadMeshVector_1.size()-1].Transposed();
    quadTemp52 = quadTemp5.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp52(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp52);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp52(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp52(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp52(i*numberOfQuadMesh+j+1, 3));
        }
    }
    quadTemp5 = D * quadMeshVector_2[quadMeshVector_2.size()-1].Transposed();
    quadTemp52 = quadTemp5.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp52(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp52);
#endif

    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            controlPoints_5.SetValue(i, j, gp_Pnt(surface_5((i - 1) * NumRows + j, 1), surface_5((i - 1) * NumRows + j, 2), surface_5((i - 1) * NumRows + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_5));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_5));

    temp56 = D * surface_5.Transposed();
    temp566 = temp56.Transposed();
    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            surface_6((i - 1) * NumRows + j, 1) = temp566((i - 1) * NumRows + j, 1);
            surface_6((i - 1) * NumRows + j, 2) = temp566((i - 1) * NumRows + j, 2);
            surface_6((i - 1) * NumRows + j, 3) = temp566((i - 1) * NumRows + j, 3) + 1;
        }
    }
    surfaceMatrixVector.push_back(surface_6);

    temp56_Second = D * surface_5_Second.Transposed();
    temp566_Second = temp56_Second.Transposed();
    for (Standard_Integer i = 1; i <= (*ctrlpts_2).ColLength(); i++) {
        for (Standard_Integer j = 1; j <= (*ctrlpts_2).RowLength(); j++) {
            surface_6_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1) = temp566_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 1);
            surface_6_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2) = temp566_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 2);
            surface_6_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) = temp566_Second((i - 1) * (*ctrlpts_2).RowLength() + j, 3) + 1;
        }
    }
    surfaceMatrixSecondVector.push_back(surface_6_Second);

#ifdef USE_QUAD_HEX_MESH
    quadTemp6 = D * quadMeshVector_1[quadMeshVector_1.size()-1].Transposed();
    quadTemp62 = quadTemp6.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp62(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_1.push_back(quadTemp62);
    quadTemp6 = D * quadMeshVector_2[quadMeshVector_2.size()-1].Transposed();
    quadTemp62 = quadTemp6.Transposed();
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadTemp62(i*numberOfQuadMesh+j+1, 3) += 1;
        }
    }
    quadMeshVector_2.push_back(quadTemp62);
    for (int i = 0; i < numberOfQuadMesh; i++) {
        for (int j = 0; j < numberOfQuadMesh; j++) {
            quadMeshPointsTotal.push_back(quadTemp62(i*numberOfQuadMesh+j+1, 1));
            quadMeshPointsTotal.push_back(quadTemp62(i*numberOfQuadMesh+j+1, 2));
            quadMeshPointsTotal.push_back(quadTemp62(i*numberOfQuadMesh+j+1, 3));
        }
    }
#endif

    for (Standard_Integer i = 1; i <= NumCols; i++) {
        for (Standard_Integer j = 1; j <= NumRows; j++) {
            controlPoints_6.SetValue(i, j, gp_Pnt(surface_6_Second((i - 1) * NumRows + j, 1), surface_6_Second((i - 1) * NumRows + j, 2), surface_6_Second((i - 1) * NumRows + j, 3)));
        }
    }
    controlPointsVector.push_back(scaleControlPoints(controlPoints_6));
    cellSurfacePointsVector.push_back(scaleControlPoints(controlPoints_6));



    for (int i = 0; i < 6; i++) {
        math_Matrix surfaceTemp(1, 3, 1, NumRows*NumCols);
        if (i % 2 == 0) {
            surfaceTemp = Ty * surfaceMatrixSecondVector[i].Transposed();
        }
        else {
            surfaceTemp = Ty * surfaceMatrixVector[i].Transposed();
        }
        math_Matrix surfaceTemp2 = surfaceTemp.Transposed();
        math_Matrix surfaceTemp3(1, controlPointsColLength * controlPointsRowLength, 1, 3);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 1) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 1) + 1;
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 2) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 2);
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 3) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 3);
            }
        }
        surfaceMatrixVector.push_back(surfaceTemp3);

#ifdef USE_QUAD_HEX_MESH
        math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
        if (i % 2 == 0) {
            quadTemp_cal1 = Ty * quadMeshVector_2[i].Transposed();
        }
        else {
            quadTemp_cal1 = Ty * quadMeshVector_1[i].Transposed();
        }
        math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 1) += 1;
            }
        }
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 1));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 2));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 3));
            }
        }
#endif

        TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                controlPointsTemp.SetValue(i, j, gp_Pnt(surfaceTemp3((i - 1) * controlPointsRowLength + j, 1), surfaceTemp3((i - 1) * controlPointsRowLength + j, 2), surfaceTemp3((i - 1) * controlPointsRowLength + j, 3)));
            }
        }
        controlPointsVector.push_back(scaleControlPoints(controlPointsTemp));
        cellSurfacePointsVector.push_back(scaleControlPoints(controlPointsTemp));
    }

    for (int i = 0; i < 6; i++) {
        math_Matrix surfaceTemp(1, 3, 1, NumRows*NumCols);
        if (i % 2 == 1) {
            surfaceTemp = Tyz * surfaceMatrixSecondVector[i].Transposed();
        }
        else {
            surfaceTemp = Tyz * surfaceMatrixVector[i].Transposed();
        }
        math_Matrix surfaceTemp2 = surfaceTemp.Transposed();
        math_Matrix surfaceTemp3(1, controlPointsColLength * controlPointsRowLength, 1, 3);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 1) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 1) + 1;
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 2) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 2) + 1;
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 3) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 3);
            }
        }
        surfaceMatrixVector.push_back(surfaceTemp3);

#ifdef USE_QUAD_HEX_MESH
        math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
        if (i % 2 == 1) {
            quadTemp_cal1 = Tyz * quadMeshVector_2[i].Transposed();
        }
        else {
            quadTemp_cal1 = Tyz * quadMeshVector_1[i].Transposed();
        }
        math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 1) += 1;
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 2) += 1;
            }
        }
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 1));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 2));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 3));
            }
        }
#endif

        TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                controlPointsTemp.SetValue(i, j, gp_Pnt(surfaceTemp3((i - 1) * controlPointsRowLength + j, 1), surfaceTemp3((i - 1) * controlPointsRowLength + j, 2), surfaceTemp3((i - 1) * controlPointsRowLength + j, 3)));
            }
        }
        controlPointsVector.push_back(scaleControlPoints(controlPointsTemp));
        cellSurfacePointsVector.push_back(scaleControlPoints(controlPointsTemp));
    }

    for (int i = 0; i < 6; i++) {
        math_Matrix surfaceTemp(1, 3, 1, NumRows*NumCols);
        if (i % 2 == 1) {
            surfaceTemp = Txy * surfaceMatrixSecondVector[i].Transposed();
        }
        else {
            surfaceTemp = Txy * surfaceMatrixVector[i].Transposed();
        }
        math_Matrix surfaceTemp2 = surfaceTemp.Transposed();
        math_Matrix surfaceTemp3(1, controlPointsColLength * controlPointsRowLength, 1, 3);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 1) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 1) + 1;
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 2) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 2);
                surfaceTemp3((i - 1) * controlPointsRowLength + j, 3) = surfaceTemp2((i - 1) * controlPointsRowLength + j, 3) - 1;
            }
        }
        surfaceMatrixVector.push_back(surfaceTemp3);

#ifdef USE_QUAD_HEX_MESH
        math_Matrix quadTemp_cal1(1, 3, 1, numberOfQuadMesh*numberOfQuadMesh);
        if (i % 2 == 1) {
            quadTemp_cal1 = Txy * quadMeshVector_2[i].Transposed();
        }
        else {
            quadTemp_cal1 = Txy * quadMeshVector_1[i].Transposed();
        }
        math_Matrix quadTemp_cal2 = quadTemp_cal1.Transposed();
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 1) += 1;
                quadTemp_cal2(i*numberOfQuadMesh+j+1, 3) -= 1;
            }
        }
        for (int i = 0; i < numberOfQuadMesh; i++) {
            for (int j = 0; j < numberOfQuadMesh; j++) {
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 1));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 2));
                quadMeshPointsTotal.push_back(quadTemp_cal2(i*numberOfQuadMesh+j+1, 3));
            }
        }
#endif

        TColgp_Array2OfPnt controlPointsTemp(1, controlPointsColLength, 1, controlPointsRowLength);
        for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
            for (Standard_Integer j = 1; j <= controlPointsRowLength; j++) {
                controlPointsTemp.SetValue(i, j, gp_Pnt(surfaceTemp3((i - 1) * controlPointsRowLength + j, 1), surfaceTemp3((i - 1) * controlPointsRowLength + j, 2), surfaceTemp3((i - 1) * controlPointsRowLength + j, 3)));
            }
        }
        controlPointsVector.push_back(scaleControlPoints(controlPointsTemp));
        cellSurfacePointsVector.push_back(controlPointsTemp);
    }

    buildCurvesGyroid2(Uknot, multiplicitiesU, controlPointsColLength, cellSurfacePointsVector, curveControlPointsVector, bSplineCurveVector);
    // make boundary surfaces
    std::vector<GeomFill_BSplineCurves> cellBoundarySurfaceVector;
    makeGapPatchesGyroid(cellBoundarySurfaceVector, bSplineCurveVector);
    // utilize the control points to build the B-Spline surface
    std::vector<Geom_BSplineSurface> basicSurfaceElementVector;
    for (int i = 0; i < 48; i++) {
        Handle(Geom_BSplineSurface) bsplineTemp = new Geom_BSplineSurface(
            controlPointsVector[i],
            *Uknot,
            *Vknot,
            (*multiplicitiesU),
            (*multiplicitiesV),
            3,
            3
        );
        cellSurfaceVector.push_back(bsplineTemp);
    }
    // duplicate and get the faces of a TPMS unit
    makeUnitFacesGyroid(basicSurfaceElementVector, cellBoundarySurfaceVector);
    // Create a BRepBuilderAPI_MakeFace to convert the B-spline surface to a face.
    std::vector<BRepBuilderAPI_MakeFace> faceMakerVector;
    std::vector<TopoDS_Face> topoFaceVector;
    int cellFaceNumber = 48;
    int internalFaceNumber = cellFaceNumber * this->getCellNumberX() * this->getCellNumberY() * this->getCellNumberZ();
    int coverNumber = 4 * 2 * (this->getCellNumberX() * this->getCellNumberY() + this->getCellNumberX() * this->getCellNumberZ() + this->getCellNumberY() * this->getCellNumberZ())
    + 2 * this->getCellNumberZ() * (2 * (this->getCellNumberX() + this->getCellNumberY()) - 2) + 2 * this->getCellNumberY() * (2 * (this->getCellNumberX() + this->getCellNumberZ()) - 2) + 2 * this->getCellNumberX() * (2 * (this->getCellNumberY() + this->getCellNumberZ()) - 2);
    int faceNumber = internalFaceNumber + coverNumber;
    std::vector<Handle(Geom_BSplineSurface)> faceVector;
    // duplicate and get the faces of the whole TPMS solid model
    makeSolidModelFacesGyroid(basicSurfaceElementVector, faceVector, cellSurfaceVector);

#ifdef USE_QUAD_HEX_MESH
    std::vector<float> quadMeshAll;
    for (int i = 0; i < this->getCellNumberX(); i++) {
        for (int j = 0; j < this->getCellNumberY(); j++) {
            for (int k = 0; k < this->getCellNumberZ(); k++) {
                for (int m = 0; m < numberOfQuadMesh*numberOfQuadMesh*48; m++) {
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3]+i*cellSizeX);
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3+1]+j*cellSizeY);
                    quadMeshAll.push_back(quadMeshPointsTotal[m*3+2]+k*cellSizeZ);
                }
            }
        }
    }
    writeHexahedralMesh(quadMeshAll, 48*this->getCellNumberX()*this->getCellNumberY()*this->getCellNumberZ(), 24, numberOfQuadMesh, "gyroid_hex.off");

    std::vector<float> quadMeshEdge;
    for (int i=0;i<faceVector.size();i++) {
        quadMeshEdge.clear();
        extractMesh(faceVector[i], quadMeshEdge, numberOfQuadMesh);
        for (int j = 0; j < numberOfQuadMesh*numberOfQuadMesh*3; j++) {
            quadMeshAll.push_back(quadMeshEdge[j]);
        }
    }
    writeMesh(quadMeshAll, 48*this->getCellNumberX()*this->getCellNumberY()*this->getCellNumberZ()+faceVector.size(), numberOfQuadMesh, "gyroid.off");
#endif

    float scaleFactor = 0.0;
    float scaleFactorx = 1.0;
    float scaleFactory = 1.0;
    float scaleFactorz = 1.0;
    for (int i = 0; i < faceNumber; i++) {
        /************    scale   **************/
        for (int j=1;j<=faceVector[i]->NbUPoles();j++) {
            for (int k=1;k<=faceVector[i]->NbVPoles();k++) {
                // faceVector[i]->SetPole(j,k,gp_Pnt((faceVector[i]->Pole(j,k)).X()*(1+scaleFactor*(faceVector[i]->Pole(j,k)).X()),
                // (faceVector[i]->Pole(j,k)).Y(), (faceVector[i]->Pole(j,k)).Z()));
                faceVector[i]->SetPole(j,k,gp_Pnt((faceVector[i]->Pole(j,k)).X()*scaleFactorx,
                (faceVector[i]->Pole(j,k)).Y()*scaleFactory, (faceVector[i]->Pole(j,k)).Z()*scaleFactorz));
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
            TopoDS_Face faceTemp = faceMakerVector[i].Face();
            topoFaceVector.push_back(faceTemp);
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
        BRepCheck_Analyzer analyzer = BRepCheck_Analyzer(*tpmsSolidModel);
        if (!analyzer.IsValid()) {
            cout << "TPMS2STEP > topo solid error" << endl;
            return false;
        }
        cout << "TPMS2STEP > topo solid valid" << endl;
        
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
        std::cerr << "TPMS2STEP > Failed to create Gyroid TPMS solid model." << std::endl;
        return false;
    }
    cout << "TPMS2STEP > Gyroid TPMS solid model generated successfully." << endl;
    return true;
}

void TPMS_G::errorControlledSampling(vector<float>& knotTemp, float* knot, int* NumRows, int* NumCols, int* numberOfPointCloud, Eigen::MatrixXf& surfacePlusMatrix,
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
        if (i == 0) {
            paramRuler.push_back((i+0.4)*1.0/(2*n-2));
        } else if (i == 2*n-3) {
            paramRuler.push_back((i+0.6)*1.0/(2*n-2));
        }
    }
    if (knot) {
        free(knot);
    }
    knotTemp.clear();
    knot = (float*)malloc((paramRuler.size()+4)*sizeof(float));
    knot[0] = 0;
    knot[1] = 0;
    knot[2] = 0;
    knot[3] = 0;
    knotTemp.push_back(0);
    knotTemp.push_back(0);
    knotTemp.push_back(0);
    knotTemp.push_back(0);
    for (int i=4;i<=paramRuler.size()-1;i++) {
        knot[i] = (i-3)*1.0/(paramRuler.size()-3);
        knotTemp.push_back((i-3)*1.0/(paramRuler.size()-3));
    }
    knot[paramRuler.size()] = 1;
    knot[paramRuler.size()+1] = 1;
    knot[paramRuler.size()+2] = 1;
    knot[paramRuler.size()+3] = 1;
    knotTemp.push_back(1);
    knotTemp.push_back(1);
    knotTemp.push_back(1);
    knotTemp.push_back(1);
    gridLength = paramRuler.size()-1;
    while (flagErrorControl) {
        iterationCount++;
        lGrid = 1.0/(paramRuler.size()-1);
        cout << "TPMS2STEP > Gyoid Error Control Iteration: " << iterationCount << endl;
        int length = paramRuler.size();// the number of control points is the same as the number of sample points
        int halfLength = (length+1)/2;
        sample2(paramRuler);
        pythonEXEC(this->modelType, this->getOffsetValue1(), this->getOffsetValue2());
        sample(paramRuler);
        pythonEXEC(this->modelType, this->getOffsetValue1(), this->getOffsetValue2());
        (*NumRows) = paramRuler.size();
        (*NumCols) = paramRuler.size();
        controlPointsPlusPIA.resize(paramRuler.size()*paramRuler.size(), 3);
        controlPointsMinusPIA.resize(paramRuler.size()*paramRuler.size(), 3);
        (*numberOfPointCloud) = paramRuler.size()*paramRuler.size();
        U.resize(*numberOfPointCloud);
        V.resize(*numberOfPointCloud);
        surfacePlusMatrix.resize((*numberOfPointCloud), 3);
        surfaceMinusMatrix.resize((*numberOfPointCloud), 3);
        // read sample points
        readSemiFundamentalPatch(surfacePlusMatrix, surfaceMinusMatrix, numberOfPointCloud, this->getOffsetValue1());
        ofs.open("../data/uv/uv_Gyroid_pia.off", ios::out);
        ofs << "OFF\n" << to_string(paramRuler.size()*paramRuler.size()) << " 0 0" << endl;
        for (int i=0;i<paramRuler.size();i++) {
            for (int j=0;j<paramRuler.size();j++) {
                U(i*paramRuler.size()+j) = j*1.0/(paramRuler.size()-1);
                V(i*paramRuler.size()+j) = i*1.0/(paramRuler.size()-1);
                ofs << j*1.0/(paramRuler.size()-1) << " " << i*1.0/(paramRuler.size()-1) << " 0" << endl;
            }
        }
        ofs.close();
        cout << "TPMS2STEP > GPU constrained-PIA" << endl;
        constrainedPIAGPU(knot, 3, (*numberOfPointCloud), (*NumRows), (*NumCols), surfacePlusMatrix, surfaceMinusMatrix, U, V,
        &controlPointsPlusPIA, &controlPointsMinusPIA, this->getOffsetValue1(), this->modelType);
        for (int i=0;i<(*NumRows)*(*NumCols);i++) {
            controlPointsPlusPIA(i,0) += 1;
            controlPointsPlusPIA(i,1) += 0.5;
            controlPointsMinusPIA(i,0) += 1;
            controlPointsMinusPIA(i,1) += 0.5;
        }
        ofs.open("controlPoint_test.txt", ios::out);
        for (int i=0;i<(*NumRows)*(*NumCols);i++) {
            ofs << controlPointsPlusPIA(i,0) << " " << controlPointsPlusPIA(i,1) << " " << controlPointsPlusPIA(i,2) << endl;
        }
        ofs.close();
        makeMinusByPlus(controlPointsPlusPIA, controlPointsMinusPIA, (*NumRows), (*NumCols));
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
        int totalRefineNumber = (paramRuler.size()-2)*2-1;
        paramRuler.clear();
        for (int i=0;i<totalRefineNumber;i++) {
            paramRuler.push_back(i*1.0/(totalRefineNumber-1));
            if (i == 0) {
                paramRuler.push_back((i+0.4)*1.0/(totalRefineNumber-1));
            } else if (i == totalRefineNumber-2) {
                paramRuler.push_back((i+0.6)*1.0/(totalRefineNumber-1));
            }
        }
        refineTimes.clear();
        refineTimes.resize(paramRuler.size()-1);
        for (int i=0;i<paramRuler.size()-1;i++) {
            refineTimes.push_back(0);
        }
        if (knot) {
            free(knot);
            knot = (float*)malloc((paramRuler.size()+4)*sizeof(float));
            knot[0] = 0;
            knot[1] = 0;
            knot[2] = 0;
            knot[3] = 0;
            for (int i=4;i<=paramRuler.size()-1;i++) {
                knot[i] = (i-3)*1.0/(paramRuler.size()-3);
            }
            knot[paramRuler.size()] = 1;
            knot[paramRuler.size()+1] = 1;
            knot[paramRuler.size()+2] = 1;
            knot[paramRuler.size()+3] = 1;
        }
    }
}

void TPMS_G::buildCurvesGyroid(TColStd_Array1OfReal* Uknot, TColStd_Array1OfInteger* multiplicitiesU, int controlPointsColLength, std::vector<TColgp_Array2OfPnt>& cellSurfacePointsVector, std::vector<TColgp_Array1OfPnt>& curveControlPointsVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector) {
    TColgp_Array1OfPnt curveCpts1(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, cellSurfacePointsVector[2].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts1);
    Handle(Geom_BSplineCurve) curveOne = new Geom_BSplineCurve(
        curveCpts1,
        (*Uknot),
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts2(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, cellSurfacePointsVector[4].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    Handle(Geom_BSplineCurve) curveTwo = new Geom_BSplineCurve(
        curveCpts2,
        (*Uknot),
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts3(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, cellSurfacePointsVector[9].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts3);
    Handle(Geom_BSplineCurve) curveThree = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts4(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, cellSurfacePointsVector[8].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    Handle(Geom_BSplineCurve) curveFour = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curveOne);
    bSplineCurveVector.push_back(curveTwo);
    bSplineCurveVector.push_back(curveThree);
    bSplineCurveVector.push_back(curveFour);

    TColgp_Array1OfPnt curveCpts5(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, cellSurfacePointsVector[6].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    Handle(Geom_BSplineCurve) curveFive = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts6(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, cellSurfacePointsVector[6].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts6);
    Handle(Geom_BSplineCurve) curveSix = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts7(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, cellSurfacePointsVector[12].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts7);
    Handle(Geom_BSplineCurve) curveSeven = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts8(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, cellSurfacePointsVector[17].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts8);
    Handle(Geom_BSplineCurve) curveEight = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts9(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts9.SetValue(i, cellSurfacePointsVector[5].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts9);
    Handle(Geom_BSplineCurve) curve9 = new Geom_BSplineCurve(
        curveCpts9,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts10(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts10.SetValue(i, cellSurfacePointsVector[5].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts10);
    Handle(Geom_BSplineCurve) curve10 = new Geom_BSplineCurve(
        curveCpts10,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts11(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts11.SetValue(i, cellSurfacePointsVector[13].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts11);
    Handle(Geom_BSplineCurve) curve11 = new Geom_BSplineCurve(
        curveCpts11,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts12(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts12.SetValue(i, cellSurfacePointsVector[14].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts12);
    Handle(Geom_BSplineCurve) curve12 = new Geom_BSplineCurve(
        curveCpts12,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts13(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts13.SetValue(i, cellSurfacePointsVector[12].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts13);
    Handle(Geom_BSplineCurve) curve13 = new Geom_BSplineCurve(
        curveCpts13,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts14(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts14.SetValue(i, cellSurfacePointsVector[10].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts14);
    Handle(Geom_BSplineCurve) curve14 = new Geom_BSplineCurve(
        curveCpts14,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts15(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts15.SetValue(i, cellSurfacePointsVector[4].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts15);
    Handle(Geom_BSplineCurve) curve15 = new Geom_BSplineCurve(
        curveCpts15,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts16(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts16.SetValue(i, cellSurfacePointsVector[3].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts16);
    Handle(Geom_BSplineCurve) curve16 = new Geom_BSplineCurve(
        curveCpts16,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts17(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts17.SetValue(i, cellSurfacePointsVector[3].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts17);
    Handle(Geom_BSplineCurve) curve17 = new Geom_BSplineCurve(
        curveCpts17,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    TColgp_Array1OfPnt curveCpts18(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts18.SetValue(i, cellSurfacePointsVector[14].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts18);
    Handle(Geom_BSplineCurve) curve18 = new Geom_BSplineCurve(
        curveCpts18,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curveFive);
    bSplineCurveVector.push_back(curveSix);
    bSplineCurveVector.push_back(curveSeven);
    bSplineCurveVector.push_back(curveEight);
    bSplineCurveVector.push_back(curve9);
    bSplineCurveVector.push_back(curve10);
    bSplineCurveVector.push_back(curve11);
    bSplineCurveVector.push_back(curve12);
    bSplineCurveVector.push_back(curve13);
    bSplineCurveVector.push_back(curve14);
    bSplineCurveVector.push_back(curve15);
    bSplineCurveVector.push_back(curve16);
    bSplineCurveVector.push_back(curve17);
    bSplineCurveVector.push_back(curve18);
}

void TPMS_G::buildCurvesGyroid2(TColStd_Array1OfReal* Uknot, TColStd_Array1OfInteger* multiplicitiesU, int controlPointsColLength, std::vector<TColgp_Array2OfPnt>& cellSurfacePointsVector, std::vector<TColgp_Array1OfPnt>& curveControlPointsVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector) {
    TColgp_Array1OfPnt curveCpts1(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts2(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts3(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts4(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts5(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts6(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts7(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts8(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts9(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts10(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts11(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts12(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts13(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts14(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts15(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts16(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts17(1, controlPointsColLength);
    TColgp_Array1OfPnt curveCpts18(1, controlPointsColLength);
    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts1.SetValue(i, cellSurfacePointsVector[26].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts1);
    Handle(Geom_BSplineCurve) curveOne = new Geom_BSplineCurve(
        curveCpts1,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts2.SetValue(i, cellSurfacePointsVector[28].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts2);
    Handle(Geom_BSplineCurve) curveTwo = new Geom_BSplineCurve(
        curveCpts2,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts3.SetValue(i, cellSurfacePointsVector[33].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts3);
    Handle(Geom_BSplineCurve) curveThree = new Geom_BSplineCurve(
        curveCpts3,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts4.SetValue(i, cellSurfacePointsVector[32].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts4);
    Handle(Geom_BSplineCurve) curveFour = new Geom_BSplineCurve(
        curveCpts4,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts5.SetValue(i, cellSurfacePointsVector[30].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts5);
    Handle(Geom_BSplineCurve) curveFive = new Geom_BSplineCurve(
        curveCpts5,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts6.SetValue(i, cellSurfacePointsVector[30].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts6);
    Handle(Geom_BSplineCurve) curveSix = new Geom_BSplineCurve(
        curveCpts6,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts7.SetValue(i, cellSurfacePointsVector[36].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts7);
    Handle(Geom_BSplineCurve) curveSeven = new Geom_BSplineCurve(
        curveCpts7,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts8.SetValue(i, cellSurfacePointsVector[41].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts8);
    Handle(Geom_BSplineCurve) curveEight = new Geom_BSplineCurve(
        curveCpts8,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts9.SetValue(i, cellSurfacePointsVector[29].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts9);
    Handle(Geom_BSplineCurve) curve9 = new Geom_BSplineCurve(
        curveCpts9,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts10.SetValue(i, cellSurfacePointsVector[29].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts10);
    Handle(Geom_BSplineCurve) curve10 = new Geom_BSplineCurve(
        curveCpts10,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts11.SetValue(i, cellSurfacePointsVector[37].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts11);
    Handle(Geom_BSplineCurve) curve11 = new Geom_BSplineCurve(
        curveCpts11,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts12.SetValue(i, cellSurfacePointsVector[38].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts12);
    Handle(Geom_BSplineCurve) curve12 = new Geom_BSplineCurve(
        curveCpts12,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts13.SetValue(i, cellSurfacePointsVector[36].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts13);
    Handle(Geom_BSplineCurve) curve13 = new Geom_BSplineCurve(
        curveCpts13,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts14.SetValue(i, cellSurfacePointsVector[34].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts14);
    Handle(Geom_BSplineCurve) curve14 = new Geom_BSplineCurve(
        curveCpts14,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts15.SetValue(i, cellSurfacePointsVector[28].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts15);
    Handle(Geom_BSplineCurve) curve15 = new Geom_BSplineCurve(
        curveCpts15,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts16.SetValue(i, cellSurfacePointsVector[27].Value(i, 1));
    }
    curveControlPointsVector.push_back(curveCpts16);
    Handle(Geom_BSplineCurve) curve16 = new Geom_BSplineCurve(
        curveCpts16,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts17.SetValue(i, cellSurfacePointsVector[27].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts17);
    Handle(Geom_BSplineCurve) curve17 = new Geom_BSplineCurve(
        curveCpts17,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    for (Standard_Integer i = 1; i <= controlPointsColLength; i++) {
        curveCpts18.SetValue(i, cellSurfacePointsVector[38].Value(1, i));
    }
    curveControlPointsVector.push_back(curveCpts18);
    Handle(Geom_BSplineCurve) curve18 = new Geom_BSplineCurve(
        curveCpts18,
        *Uknot,
        (*multiplicitiesU),
        3
    );

    bSplineCurveVector.push_back(curveOne);
    bSplineCurveVector.push_back(curveTwo);
    bSplineCurveVector.push_back(curveThree);
    bSplineCurveVector.push_back(curveFour);
    bSplineCurveVector.push_back(curveFive);
    bSplineCurveVector.push_back(curveSix);
    bSplineCurveVector.push_back(curveSeven);
    bSplineCurveVector.push_back(curveEight);
    bSplineCurveVector.push_back(curve9);
    bSplineCurveVector.push_back(curve10);
    bSplineCurveVector.push_back(curve11);
    bSplineCurveVector.push_back(curve12);
    bSplineCurveVector.push_back(curve13);
    bSplineCurveVector.push_back(curve14);
    bSplineCurveVector.push_back(curve15);
    bSplineCurveVector.push_back(curve16);
    bSplineCurveVector.push_back(curve17);
    bSplineCurveVector.push_back(curve18);
}

void TPMS_G::readSemiFundamentalPatch(Eigen::MatrixXf& surfacePlusMatrix, Eigen::MatrixXf& surfaceMinusMatrix, int* numberOfPointCloud, float offsetValue) {
    string tempString;
    ifstream ifs;
    ifs = ifstream("../data/SFP/Gyroid/SFP_Gyroid_offset_" + to_string(int(offsetValue * (-100) + 1 - 1)) + ".off", ios::in);
    ifs >> tempString;
    int vertex, edge, face;
    float a, b, c;
    ifs >> vertex >> edge >> face;
    (*numberOfPointCloud) = vertex;
    for (Standard_Integer i = 0; i < vertex; i++) {
        ifs >> a >> b >> c;
        surfacePlusMatrix(i, 0) = a;
        surfacePlusMatrix(i, 1) = b;
        surfacePlusMatrix(i, 2) = c;
    }
    ifs.close();
    ifs = ifstream("../data/SFP/Gyroid/SFP_Gyroid_offset_" + to_string(int((offsetValue * (-100) + 1 - 1) * (-1))) + ".off", ios::in);
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

void TPMS_G::makeGapPatchesGyroid(std::vector<GeomFill_BSplineCurves>& cellBoundarySurfaceVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector) {
    int edgeBase = 18;
    // x axis
    gp_Pnt startPointFirst = bSplineCurveVector[4]->StartPoint();
    gp_Pnt startPointSecond = bSplineCurveVector[4 + edgeBase]->StartPoint();
    gp_Pnt endPointFirst = bSplineCurveVector[4]->EndPoint();
    gp_Pnt endPointSecond = bSplineCurveVector[4 + edgeBase]->EndPoint();

    TColgp_Array1OfPnt curveControl(1, 2);
    curveControl.SetValue(1, startPointFirst);
    curveControl.SetValue(2, startPointSecond);
    TColStd_Array1OfReal curveKnot(1, 2);
    curveKnot.SetValue(1, 0);
    curveKnot.SetValue(2, 1);
    TColStd_Array1OfInteger multiplicitiesCurve(1, 2);
    multiplicitiesCurve.SetValue(1, 2);
    multiplicitiesCurve.SetValue(2, 2);

    Handle(Geom_BSplineCurve) curveFirst = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );
    curveControl.SetValue(1, endPointFirst);
    curveControl.SetValue(2, endPointSecond);
    Handle(Geom_BSplineCurve) curveSecond = new Geom_BSplineCurve(
        curveControl,
        curveKnot,
        multiplicitiesCurve,
        1
    );

    GeomFill_BSplineCurves aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[4], curveFirst, bSplineCurveVector[4 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);

    startPointFirst = bSplineCurveVector[5]->StartPoint();
    startPointSecond = bSplineCurveVector[5 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[5]->EndPoint();
    endPointSecond = bSplineCurveVector[5 + edgeBase]->EndPoint();

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
    GeomFill_BSplineCurves aFillSurface2 = GeomFill_BSplineCurves(bSplineCurveVector[5], curveFirst, bSplineCurveVector[5 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface2);

    startPointFirst = bSplineCurveVector[6]->StartPoint();
    startPointSecond = bSplineCurveVector[6 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[6]->EndPoint();
    endPointSecond = bSplineCurveVector[6 + edgeBase]->EndPoint();

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

    GeomFill_BSplineCurves aFillSurface3 = GeomFill_BSplineCurves(bSplineCurveVector[6], curveFirst, bSplineCurveVector[6 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface3);

    startPointFirst = bSplineCurveVector[7]->StartPoint();
    startPointSecond = bSplineCurveVector[7 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[7]->EndPoint();
    endPointSecond = bSplineCurveVector[7 + edgeBase]->EndPoint();

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
    GeomFill_BSplineCurves aFillSurface4 = GeomFill_BSplineCurves(bSplineCurveVector[7], curveFirst, bSplineCurveVector[7 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface4);
    // y axis
    startPointFirst = bSplineCurveVector[8]->StartPoint();
    startPointSecond = bSplineCurveVector[8 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[8]->EndPoint();
    endPointSecond = bSplineCurveVector[8 + edgeBase]->EndPoint();

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
    aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[8], curveFirst, bSplineCurveVector[8 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);

    startPointFirst = bSplineCurveVector[9]->StartPoint();
    startPointSecond = bSplineCurveVector[9 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[9]->EndPoint();
    endPointSecond = bSplineCurveVector[9 + edgeBase]->EndPoint();

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
    aFillSurface2 = GeomFill_BSplineCurves(bSplineCurveVector[9], curveFirst, bSplineCurveVector[9 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface2);

    startPointFirst = bSplineCurveVector[10]->StartPoint();
    startPointSecond = bSplineCurveVector[10 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[10]->EndPoint();
    endPointSecond = bSplineCurveVector[10 + edgeBase]->EndPoint();

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
    aFillSurface3 = GeomFill_BSplineCurves(bSplineCurveVector[10], curveFirst, bSplineCurveVector[10 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface3);

    startPointFirst = bSplineCurveVector[11]->StartPoint();
    startPointSecond = bSplineCurveVector[11 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[11]->EndPoint();
    endPointSecond = bSplineCurveVector[11 + edgeBase]->EndPoint();

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
    aFillSurface4 = GeomFill_BSplineCurves(bSplineCurveVector[11], curveFirst, bSplineCurveVector[11 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface4);

    // z axis
    startPointFirst = bSplineCurveVector[0]->StartPoint();
    startPointSecond = bSplineCurveVector[edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[0]->EndPoint();
    endPointSecond = bSplineCurveVector[edgeBase]->EndPoint();

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
    aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[0], curveFirst, bSplineCurveVector[edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);

    startPointFirst = bSplineCurveVector[1]->StartPoint();
    startPointSecond = bSplineCurveVector[1 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[1]->EndPoint();
    endPointSecond = bSplineCurveVector[1 + edgeBase]->EndPoint();

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
    aFillSurface2 = GeomFill_BSplineCurves(bSplineCurveVector[1], curveFirst, bSplineCurveVector[1 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface2);

    startPointFirst = bSplineCurveVector[2]->StartPoint();
    startPointSecond = bSplineCurveVector[2 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[2]->EndPoint();
    endPointSecond = bSplineCurveVector[2 + edgeBase]->EndPoint();

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
    aFillSurface3 = GeomFill_BSplineCurves(bSplineCurveVector[2], curveFirst, bSplineCurveVector[2 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface3);

    startPointFirst = bSplineCurveVector[3]->StartPoint();
    startPointSecond = bSplineCurveVector[3 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[3]->EndPoint();
    endPointSecond = bSplineCurveVector[3 + edgeBase]->EndPoint();

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
    aFillSurface4 = GeomFill_BSplineCurves(bSplineCurveVector[3], curveFirst, bSplineCurveVector[3 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface4);

    // x y axis
    startPointFirst = bSplineCurveVector[12]->StartPoint();
    startPointSecond = bSplineCurveVector[12 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[12]->EndPoint();
    endPointSecond = bSplineCurveVector[12 + edgeBase]->EndPoint();

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
    aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[12], curveFirst, bSplineCurveVector[12 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);

    // x z axis
    startPointFirst = bSplineCurveVector[13]->StartPoint();
    startPointSecond = bSplineCurveVector[13 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[13]->EndPoint();
    endPointSecond = bSplineCurveVector[13 + edgeBase]->EndPoint();

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
    aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[13], curveFirst, bSplineCurveVector[13 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);

    // y z axis
    startPointFirst = bSplineCurveVector[14]->StartPoint();
    startPointSecond = bSplineCurveVector[14 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[14]->EndPoint();
    endPointSecond = bSplineCurveVector[14 + edgeBase]->EndPoint();

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
    aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[14], curveFirst, bSplineCurveVector[14 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);

    // x -z axis
    startPointFirst = bSplineCurveVector[15]->StartPoint();
    startPointSecond = bSplineCurveVector[15 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[15]->EndPoint();
    endPointSecond = bSplineCurveVector[15 + edgeBase]->EndPoint();

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
    aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[15], curveFirst, bSplineCurveVector[15 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);

    // x -y axis
    startPointFirst = bSplineCurveVector[16]->StartPoint();
    startPointSecond = bSplineCurveVector[16 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[16]->EndPoint();
    endPointSecond = bSplineCurveVector[16 + edgeBase]->EndPoint();

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
    aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[16], curveFirst, bSplineCurveVector[16 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);

    // y -z axis
    startPointFirst = bSplineCurveVector[17]->StartPoint();
    startPointSecond = bSplineCurveVector[17 + edgeBase]->StartPoint();
    endPointFirst = bSplineCurveVector[17]->EndPoint();
    endPointSecond = bSplineCurveVector[17 + edgeBase]->EndPoint();

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
    aFillSurface = GeomFill_BSplineCurves(bSplineCurveVector[17], curveFirst, bSplineCurveVector[17 + edgeBase], curveSecond, GeomFill_StretchStyle);
    cellBoundarySurfaceVector.push_back(aFillSurface);
}

void TPMS_G::makeUnitFacesGyroid(std::vector<Geom_BSplineSurface>& basicSurfaceElementVector, std::vector<GeomFill_BSplineCurves>& cellBoundarySurfaceVector) {
    for (int i = 0; i < 4; i++) {
        Handle(Geom_BSplineSurface) surf = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i].Surface()->Translated(gp_Vec((-1)*cellSizeX, 0, 0)));
        basicSurfaceElementVector.push_back(*surf);
    }
    for (int i = 4; i < 8; i++) {
        Handle(Geom_BSplineSurface) surf = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i].Surface()->Translated(gp_Vec(0, (-1)*cellSizeY, 0)));
        basicSurfaceElementVector.push_back(*surf);
    }
    for (int i = 8; i < 12; i++) {
        Handle(Geom_BSplineSurface) surf = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[i].Surface()->Translated(gp_Vec(0, 0, (-1)*cellSizeZ)));
        basicSurfaceElementVector.push_back(*surf);
    }
    Handle(Geom_BSplineSurface) surfTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[12].Surface()->Translated(gp_Vec((-1)*cellSizeX, (-1)*cellSizeY, 0)));
    basicSurfaceElementVector.push_back(*surfTemp);
    surfTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[13].Surface()->Translated(gp_Vec((-1)*cellSizeX, 0, (-1)*cellSizeZ)));
    basicSurfaceElementVector.push_back(*surfTemp);
    surfTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[14].Surface()->Translated(gp_Vec(0, (-1)*cellSizeY, (-1)*cellSizeZ)));
    basicSurfaceElementVector.push_back(*surfTemp);
    surfTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[15].Surface()->Translated(gp_Vec(cellSizeX, 0, (-1)*cellSizeZ)));
    basicSurfaceElementVector.push_back(*surfTemp);
    surfTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[16].Surface()->Translated(gp_Vec(cellSizeX, (-1)*cellSizeY, 0)));
    basicSurfaceElementVector.push_back(*surfTemp);
    surfTemp = Handle(Geom_BSplineSurface)::DownCast(cellBoundarySurfaceVector[17].Surface()->Translated(gp_Vec(0, 0, 0)));
    basicSurfaceElementVector.push_back(*surfTemp);
}

void TPMS_G::makeSolidModelFacesGyroid(std::vector<Geom_BSplineSurface>& basicSurfaceElementVector, std::vector<Handle(Geom_BSplineSurface)>& faceVector, std::vector<Handle(Geom_BSplineSurface)>& cellSurfaceVector) {
    for (int i = 0; i < 4; i++) {
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

    for (int i = 4; i < 8; i++) {
        Handle(Geom_BSplineSurface) tempFace = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[i].Copy());
        for (int i = 0; i < this->getCellNumberX(); i++) {
            for (int j = 0; j < this->getCellNumberZ(); j++) {
                Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace->Translated(gp_Vec(i * cellSizeX, 0, j * cellSizeZ)));
                faceVector.push_back(tempTranslatedFace);
                Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, this->getCellNumberY() * cellSizeY, 0)));
                faceVector.push_back(tempTranslatedFace2);
            }
        }
    }

    for (int i = 8; i < 12; i++) {
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

    Handle(Geom_BSplineSurface) tempFace12 = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[12].Copy());
    for (int i = 0; i < this->getCellNumberZ(); i++) {
        for (int j = 0; j < this->getCellNumberX(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace12->Translated(gp_Vec(j * cellSizeX, 0, i * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(cellSizeX, this->getCellNumberY() * cellSizeY, 0)));
            faceVector.push_back(tempTranslatedFace2);
        }
        for (int j = 1; j < this->getCellNumberY(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace12->Translated(gp_Vec(0, j * cellSizeY, i * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(this->getCellNumberX() * cellSizeX, 0, 0)));
            faceVector.push_back(tempTranslatedFace2);
        }
    }

    Handle(Geom_BSplineSurface) tempFace13 = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[13].Copy());
    for (int i = 0; i < this->getCellNumberY(); i++) {
        for (int j = 0; j < this->getCellNumberX(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace13->Translated(gp_Vec(j * cellSizeX, i * cellSizeY, 0)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(cellSizeX, 0, this->getCellNumberZ() * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace2);
        }
        for (int j = 1; j < this->getCellNumberZ(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace13->Translated(gp_Vec(0, i * cellSizeY, j * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(this->getCellNumberX() * cellSizeX, 0, 0)));
            faceVector.push_back(tempTranslatedFace2);
        }
    }

    Handle(Geom_BSplineSurface) tempFace14 = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[14].Copy());
    for (int i = 0; i < this->getCellNumberX(); i++) {
        for (int j = 0; j < this->getCellNumberY(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace14->Translated(gp_Vec(i * cellSizeX, j * cellSizeY, 0)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, cellSizeY, this->getCellNumberZ() * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace2);
        }
        for (int j = 1; j < this->getCellNumberZ(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace14->Translated(gp_Vec(i * cellSizeX, 0, j * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, this->getCellNumberY() * cellSizeY, 0)));
            faceVector.push_back(tempTranslatedFace2);
        }
    }

    Handle(Geom_BSplineSurface) tempFace15 = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[15].Copy());
    for (int i = 0; i < this->getCellNumberY(); i++) {
        for (int j = 0; j < this->getCellNumberX(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace15->Translated(gp_Vec(j * cellSizeX, i * cellSizeY, 0)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec((-1) * cellSizeX, 0, this->getCellNumberZ() * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace2);
        }
        for (int j = 1; j < this->getCellNumberZ(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace15->Translated(gp_Vec((-1) * cellSizeX, i * cellSizeY, j * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(this->getCellNumberX() * cellSizeX, 0, 0)));
            faceVector.push_back(tempTranslatedFace2);
        }
    }

    Handle(Geom_BSplineSurface) tempFace16 = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[16].Copy());
    for (int i = 0; i < this->getCellNumberZ(); i++) {
        for (int j = 0; j < this->getCellNumberX(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace16->Translated(gp_Vec(j * cellSizeX, 0, i * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec((-1) * cellSizeX, this->getCellNumberY() * cellSizeY, 0)));
            faceVector.push_back(tempTranslatedFace2);
        }
        for (int j = 1; j < this->getCellNumberY(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace16->Translated(gp_Vec((-1) * cellSizeX, j * cellSizeY, i * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(this->getCellNumberX() * cellSizeX, 0, 0)));
            faceVector.push_back(tempTranslatedFace2);
        }
    }

    Handle(Geom_BSplineSurface) tempFace17 = Handle(Geom_BSplineSurface)::DownCast(basicSurfaceElementVector[17].Copy());
    for (int i = 0; i < this->getCellNumberX(); i++) {
        for (int j = 0; j < this->getCellNumberY(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace17->Translated(gp_Vec(i * cellSizeX, j * cellSizeY, 0)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, (-1) * cellSizeY, this->getCellNumberZ() * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace2);
        }
        for (int j = 1; j < this->getCellNumberZ(); j++) {
            Handle(Geom_BSplineSurface) tempTranslatedFace = Handle(Geom_BSplineSurface)::DownCast(tempFace17->Translated(gp_Vec(i * cellSizeX, (-1) * cellSizeY, j * cellSizeZ)));
            faceVector.push_back(tempTranslatedFace);
            Handle(Geom_BSplineSurface) tempTranslatedFace2 = Handle(Geom_BSplineSurface)::DownCast(tempTranslatedFace->Translated(gp_Vec(0, this->getCellNumberY() * cellSizeY, 0)));
            faceVector.push_back(tempTranslatedFace2);
        }
    }

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
}

void TPMS_G::makeMinusByPlus(Eigen::MatrixXf& controlPointsPlusPIA, Eigen::MatrixXf& controlPointsMinusPIA, int NumRows, int NumCols) {
    // use T_{3g} as the rigid transformation matrix
    Eigen::MatrixXf temp(NumRows*NumCols, 3);
    for (int i=0;i<NumRows*NumCols;i++) {
        temp(i,0) = 0 * controlPointsPlusPIA(i,0) + (-1) * controlPointsPlusPIA(i,1) + 0 * controlPointsPlusPIA(i,2) + 1.5;
        temp(i,1) = 1 * controlPointsPlusPIA(i,0) + 0 * controlPointsPlusPIA(i,1) + 0 * controlPointsPlusPIA(i,2) - 0.5;
        temp(i,2) = 0 * controlPointsPlusPIA(i,0) + 0 * controlPointsPlusPIA(i,1) + (-1) * controlPointsPlusPIA(i,2) + 0.5;
    }
    // re-order
    for (int i=0;i<NumRows;i++) {
        for (int j=0;j<NumCols;j++) {
            controlPointsMinusPIA(i*NumCols+j,0) = temp(j*NumCols+(NumRows-1-i));
        }
    }
}