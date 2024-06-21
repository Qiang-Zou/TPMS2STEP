#ifndef TPMS_P_H
#define TPMS_P_H

#include "sample.h"
#include "pythonEXEC.h"
#include "TPMS_Solid.h"

using namespace std;

class TPMS_P : public TPMS_Solid {
    private:
        
    public:
        TPMS_P();
        ~TPMS_P();

        const string modelType = "SchwarzP";
        const string outputSTEPFileName = "SchwarzP.STEP";
        int minControlPointNumberPerEdge;


        virtual bool constructSolidModel();

        virtual void errorControlledSampling(vector<float>& knotTemp, float* knot, int* NumRows, int* NumCols, int* numberOfPointCloud, Eigen::MatrixXf& surfacePlusMatrix,
        Eigen::MatrixXf& surfaceMinusMatrix, Eigen::VectorXf& U, Eigen::VectorXf& V, Eigen::MatrixXf& controlPointsPlusPIA, Eigen::MatrixXf& controlPointsMinusPIA,
        int n, float epsilon_2, std::vector<float>& paramRuler);

        void readSemiFundamentalPatch(Eigen::MatrixXf& surfacePlusMatrix, Eigen::MatrixXf& surfaceMinusMatrix, int* numberOfPointCloud, float offsetValue2);

        void buildCurvesSchwarzP(std::vector<Handle(Geom_BSplineCurve)>& newCurve, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector, std::vector<Geom_BSplineSurface>& basicSurfaceElementVector);

};

#endif