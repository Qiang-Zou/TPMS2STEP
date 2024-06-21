// QSurfaceMesh.h: interface for the QMeshPatch class.
//
//////////////////////////////////////////////////////////////////////

#ifndef QSurfaceMesh_h
#define QSurfaceMesh_h

// undef OpenMesh does not compile with min/max macros active! Please add a define NOMINMAX to your compiler flags or add #undef max before including OpenMesh headers !")
// error max macro active
#undef min
#undef max

#include <vector>
#include <algorithm>
#include <memory>

#include "Geometry.h"
#include "QMesh/QMeshPatch.h"
#include "OpenMesh.h"

class QSurfaceMesh : public QMeshPatch, public QOpenMesh, public Geometry {
public:
    QSurfaceMesh() {QSurfaceMesh::QMeshPatch();}
    ~QSurfaceMesh() {}

    // sync data
    bool meshPatch2OpenMesh();
    bool openMesh2MeshPatch();
    void deleteOpenMeshPart();

    GeometryTypes getGeometryType() {return SurfaceMeshQMeshPatch;}
    std::shared_ptr<Geometry> getEntities2Slice(bool isIncremental, double height);
};

#endif
