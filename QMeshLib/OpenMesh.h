//
//////////////////////////////////////////////////////////////////////

#ifndef OPENMESH_H
#define OPENMESH_H

#include <vector>
#include <algorithm>
using namespace std;

#include "../GLKLib/GLKObList.h"
#include "Geometry.h"

#include "OpenMesh/Core/Mesh/Traits.hh"
#include "OpenMesh/Core/Mesh/PolyMeshT.hh"
#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"
#include "OpenMesh/Core/Geometry/EigenVectorT.hh"

struct MeshTraits : OpenMesh::DefaultTraits {
    using Scalar = double;
    using Point = Eigen::Vector3d;
    using Normal = Eigen::Vector3d;
    using TexCoord2D = Eigen::Vector2d;
    //enable standart properties
    VertexAttributes (OpenMesh::Attributes::Status|OpenMesh::Attributes::Normal|OpenMesh::Attributes::Color);
    HalfedgeAttributes (OpenMesh::Attributes::Status|OpenMesh::Attributes::PrevHalfedge);
    FaceAttributes (OpenMesh::Attributes::Status|OpenMesh::Attributes::Normal|OpenMesh::Attributes::Color);
    EdgeAttributes (OpenMesh::Attributes::Status|OpenMesh::Attributes::Color);
};

//class QOpenMesh : public OpenMesh::PolyMesh_ArrayKernelT<OpenMesh::DefaultTraits>,
//                 public Geometry, public GLKObList {
//public:
//    QOpenMesh() {}
//    ~QOpenMesh();

//    GeometryTypes GetGeometryType() {return SurfaceMeshOpenMesh;}
//    std::shared_ptr<Geometry> GetEntities2Slice(bool isIncremental, double height);

//};

class QOpenMesh : public OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> {
public:
    GeometryTypes GetGeometryType();

};

#endif
