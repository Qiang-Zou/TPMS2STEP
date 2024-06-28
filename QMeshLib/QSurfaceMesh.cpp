#include "QSurfaceMesh.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "QMesh/QMeshPatch.h"
#include "QMesh/QMeshFace.h"
#include "QMesh/QMeshEdge.h"
#include "QMesh/QMeshNode.h"

#include <OpenMesh/Core/IO/MeshIO.hh>


std::shared_ptr<Geometry> QSurfaceMesh::getEntities2Slice(bool isIncremental, double height)
{
    return nullptr;
}

bool QSurfaceMesh::meshPatch2OpenMesh()
{
    // write the mesh to a temporary place
    std::string path = "./out.obj";
    this->outputOBJFile(&path[0]);

    // load it to the OpenMesh
    if (!OpenMesh::IO::read_mesh(*this, "./out.obj")) {
        std::cerr << "QSurfaceMesh > error reading file" << endl;
        return false;
    }

    // delete the temporary file
    if(remove("./out.obj") != 0) {
        std::cerr << "QSurfaceMesh > Error deleting file" << endl;
        return false;
    }

    return true;
}

bool QSurfaceMesh::openMesh2MeshPatch()
{
    // write the mesh to a temporary place
    if (!OpenMesh::IO::write_mesh(*this, "./out.obj")) {
        std::cerr << "QSurfaceMesh > error writing file" << endl;
        return false;
    }

    // load it to the QMeshPatch
    std::string path = "./out.obj";
    if (!this->inputOBJFile(&path[0])) {
        std::cerr << "QSurfaceMesh > error reading file" << endl;
        return false;
    }

    // delete the temporary file
    if(remove("./out.obj") != 0) {
        std::cerr << "QSurfaceMesh > Error deleting file" << endl;
        return false;
    }

    return true;
}

void QSurfaceMesh::deleteOpenMeshPart()
{
    for (auto iter(this->vertices_begin()); iter != this->vertices_end(); ++iter) {
        this->delete_vertex(iter.handle(), true);
    }
    this->garbage_collection();
}
