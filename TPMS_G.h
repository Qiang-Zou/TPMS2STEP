#ifndef TPMS_G_H
#define TPMS_G_H

#include "sample.h"
#include "pythonEXEC.h"
#include "TPMS_Solid.h"

class TPMS_G : public TPMS_Solid {
    private:
        
    public:
        TPMS_G();
        ~TPMS_G();

        const string modelType = "Gyroid";
        const string outputSTEPFileName = "Gyroid.STEP";
        int minControlPointNumberPerEdge;

        virtual bool constructSolidModel();

        virtual void errorControlledSampling(vector<float>& knotTemp, float* knot, int* NumRows, int* NumCols, int* numberOfPointCloud, Eigen::MatrixXf& surfacePlusMatrix,
        Eigen::MatrixXf& surfaceMinusMatrix, Eigen::VectorXf& U, Eigen::VectorXf& V, Eigen::MatrixXf& controlPointsPlusPIA, Eigen::MatrixXf& controlPointsMinusPIA,
        int n, float epsilon_2, std::vector<float>& paramRuler);

        void buildCurvesGyroid(TColStd_Array1OfReal* Uknot, TColStd_Array1OfInteger* multiplicitiesU, int controlPointsColLength, std::vector<TColgp_Array2OfPnt>& cellSurfacePointsVector,
        std::vector<TColgp_Array1OfPnt>& curveControlPointsVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector);
        
        void buildCurvesGyroid2(TColStd_Array1OfReal* Uknot, TColStd_Array1OfInteger* multiplicitiesU, int controlPointsColLength, std::vector<TColgp_Array2OfPnt>& cellSurfacePointsVector,
        std::vector<TColgp_Array1OfPnt>& curveControlPointsVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector);
        
        void readSemiFundamentalPatch(Eigen::MatrixXf& surfacePlusMatrix, Eigen::MatrixXf& surfaceMinusMatrix, int* numberOfPointCloud, float offsetValue2);

        void makeGapPatchesGyroid(std::vector<GeomFill_BSplineCurves>& cellBoundarySurfaceVector, std::vector<Handle(Geom_BSplineCurve)>& bSplineCurveVector);

        void makeUnitFacesGyroid(std::vector<Geom_BSplineSurface>& basicSurfaceElementVector, std::vector<GeomFill_BSplineCurves>& cellBoundarySurfaceVector);

        void makeSolidModelFacesGyroid(std::vector<Geom_BSplineSurface>& basicSurfaceElementVector, std::vector<Handle(Geom_BSplineSurface)>& faceVector, std::vector<Handle(Geom_BSplineSurface)>& cellSurfaceVector);

        void makeMinusByPlus(Eigen::MatrixXf &controlPointsPlusPIA, Eigen::MatrixXf &controlPointsMinusPIA, int NumRows, int NumCols);
};

#endif