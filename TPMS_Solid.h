#ifndef TPMS_SOLID_H
#define TPMS_SOLID_H

#include <iostream>
#include <fstream>
#include <string>

// Eigen
#include <Eigen/Core>

// OCCT
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <GeomFill_BSplineCurves.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>
#include <STEPControl_StepModelType.hxx>
#include <STEPControl_Writer.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Builder.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopTools_HSequenceOfShape.hxx>

// local head files
#include "piaCuda.h"

using namespace std;

class TPMS_Solid
{
private:
    float tolerance;
    float offsetValue1;
    float offsetValue2;
    int cellNumberX;
    int cellNumberY;
    int cellNumberZ;
    
public:
    TPMS_Solid();
    ~TPMS_Solid();
    TopoDS_Solid* tpmsSolidModel;

    // argument for error control
    int maxIteration = 3;
    double cellScaleX, cellScaleY, cellScaleZ;
    double cellSizeX, cellSizeY, cellSizeZ;

    void setTolerance(float);
    float getTolerance();

    void setOffsetValue(float, float);
    float getOffsetValue1();
    float getOffsetValue2();

    void setCellNumbers(int, int, int);
    int getCellNumberX();
    int getCellNumberY();
    int getCellNumberZ();

    virtual bool constructSolidModel(){return false;};

    virtual void errorControlledSampling(vector<float>& knotTemp, float* knot, int* NumRows, int* NumCols, int* numberOfPointCloud, Eigen::MatrixXf& surfacePlusMatrix,
    Eigen::MatrixXf& surfaceMinusMatrix, Eigen::VectorXf& U, Eigen::VectorXf& V, Eigen::MatrixXf& controlPointsPlusPIA, Eigen::MatrixXf& controlPointsMinusPIA,
    int n, float epsilon_2, std::vector<float>& paramRuler){};

    // write to step files
    int writeSTEP(string outputFileName, TopoDS_Solid* solid);
    int writeSTEP(string outputFileName, int faceNumber, std::vector<TopoDS_Face>& topoFaces);
    
    // auxillary functions
    math_Matrix Rz(float);
    math_Matrix Ry(float);
    math_Matrix Rx(float);
    TColgp_Array2OfPnt scaleControlPoints(TColgp_Array2OfPnt ctrlpts);
};

extern void extractMesh(Handle_Geom_BSplineSurface, std::vector<float>&, int);
extern void writeMesh(std::vector<float>&, int, int, std::string);
extern void writeHexahedralMesh(std::vector<float>&, int, int, int, std::string);

#endif