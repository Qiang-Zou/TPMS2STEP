#include "TPMS_Solid.h"

TPMS_Solid::TPMS_Solid()
{
    tpmsSolidModel = NULL;
}

TPMS_Solid::~TPMS_Solid()
{
}

void TPMS_Solid::setTolerance(float tol) {
    this->tolerance = tol;
}

float TPMS_Solid::getTolerance() {
    return this->tolerance;
}

void TPMS_Solid::setOffsetValue(float ofv1, float ofv2) {
    this->offsetValue1 = ofv1;
    this->offsetValue2 = ofv2;
}

float TPMS_Solid::getOffsetValue1() {
    return this->offsetValue1;
}

float TPMS_Solid::getOffsetValue2() {
    return this->offsetValue2;
}

void TPMS_Solid::setCellNumbers(int x, int y, int z) {
    this->cellNumberX = x;
    this->cellNumberY = y;
    this->cellNumberZ = z;
}

int TPMS_Solid::getCellNumberX() {
    return this->cellNumberX;
}

int TPMS_Solid::getCellNumberY(){
    return this->cellNumberY;
}

int TPMS_Solid::getCellNumberZ(){
    return this->cellNumberZ;
}

TColgp_Array2OfPnt TPMS_Solid::scaleControlPoints(TColgp_Array2OfPnt ctrlpts) {
    for (int i=1;i<=ctrlpts.ColLength();i++) {
        for (int j=1;j<=ctrlpts.RowLength();j++) {
            float x = ctrlpts.Value(i,j).X()*cellScaleX;
            float y = ctrlpts.Value(i,j).Y()*cellScaleY;
            float z = ctrlpts.Value(i,j).Z()*cellScaleZ;
            ctrlpts.SetValue(i, j, gp_Pnt(x,y,z));
        }
    }
    return ctrlpts;
}

math_Matrix TPMS_Solid::Rx(float phi) {
    math_Matrix data(1, 3, 1, 3, 0);
    data(1, 1) = 1;
    data(2, 2) = cos(phi); data(2, 3) = (-1) * sin(phi);
    data(3, 2) = sin(phi); data(3, 3) = cos(phi);
    return data;
}

math_Matrix TPMS_Solid::Ry(float phi) {
    math_Matrix data(1, 3, 1, 3, 0);
    data(1, 1) = cos(phi); data(1, 3) = sin(phi);
    data(2, 2) = 1;
    data(3, 1) = (-1) * sin(phi); data(3, 3) = cos(phi);
    return data;
}

math_Matrix TPMS_Solid::Rz(float phi) {
    math_Matrix data(1, 3, 1, 3, 0);
    data(1, 1) = cos(phi); data(1, 2) = (-1) * sin(phi);
    data(2, 1) = sin(phi); data(2, 2) = cos(phi);
    data(3, 3) = 1;
    return data;
}

int TPMS_Solid::writeSTEP(string outputFileName, int faceNumber, std::vector<TopoDS_Face>& topoFaces) {
    Handle(TopTools_HSequenceOfShape) aHSequenceOfShape = new TopTools_HSequenceOfShape;
    aHSequenceOfShape->Clear();
    for (int i = 0; i < faceNumber; i++) {
        aHSequenceOfShape->Append(topoFaces[i]);
    }
    Standard_CString fileName = outputFileName.c_str();
    // save file
    IFSelect_ReturnStatus status;
    STEPControl_StepModelType aValue = STEPControl_AsIs;
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
    return 0;
}

int TPMS_Solid::writeSTEP(string outputFileName, TopoDS_Solid* solid) {
    Handle(TopTools_HSequenceOfShape) aHSequenceOfShape = new TopTools_HSequenceOfShape;
    aHSequenceOfShape->Clear();
    aHSequenceOfShape->Append(*solid);
    Standard_CString fileName = outputFileName.c_str();
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
    return 0;
}


