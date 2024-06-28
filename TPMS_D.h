#ifndef TPMS_D_H
#define TPMS_D_H

#include "sample.h"
#include "pythonEXEC.h"
#include "TPMS_Solid.h"

using namespace std;

class TPMS_D : public TPMS_Solid {
    private:
        
    public:
        TPMS_D();
        ~TPMS_D();

        const string modelType = "Diamond";
        const string outputSTEPFileName = "Diamond.STEP";
        int minControlPointNumberPerEdge;

        virtual bool constructSolidModel();

        void errorControlledSampling(vector<float>& knotTemp, float* knot, int* NumRows, int* NumCols, int* numberOfPointCloud, Eigen::MatrixXf& surfacePlusMatrix,
        Eigen::MatrixXf& surfaceMinusMatrix, Eigen::VectorXf& U, Eigen::VectorXf& V, Eigen::MatrixXf& controlPointsPlusPIA, Eigen::MatrixXf& controlPointsMinusPIA,
        int n, float epsilon_2, std::vector<float>& paramRuler);

        void buildCurvesDiamond(TColStd_Array1OfReal* Uknot, TColStd_Array1OfInteger* multiplicitiesU, int controlPointsColLength, std::vector<TColgp_Array2OfPnt>& cellSurfacePointsVector,
        std::vector<TColgp_Array1OfPnt>& curveControlPointsVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector, std::vector<TColgp_Array2OfPnt>& controlPointsTotalVector);
        
        void buildCurvesDiamond2(TColStd_Array1OfReal* Uknot, TColStd_Array1OfInteger* multiplicitiesU, int controlPointsColLength, std::vector<TColgp_Array2OfPnt>& cellSurfacePointsVector,
        std::vector<TColgp_Array1OfPnt>& curveControlPointsVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector, std::vector<TColgp_Array2OfPnt>& controlPointsTotalVector);
        
        void readSemiFundamentalPatch(Eigen::MatrixXf& surfacePlusMatrix, Eigen::MatrixXf& surfaceMinusMatrix, int* numberOfPointCloud, float offsetValue2);

        void makeSolidModelFacesDiamond(std::vector<Handle(Geom_BSplineSurface)>& faceVector, std::vector<Handle(Geom_BSplineSurface)>& cellSurfaceVector, std::vector<Handle(Geom_BSplineSurface)>& cellBoundarySurfaceVector);
};

#endif