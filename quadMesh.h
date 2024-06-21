#ifndef QUAD_MESH_H
#define QUAD_MESH_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>

#include <Geom_BSplineSurface.hxx>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>
#include <typeinfo>
#include <cstring>

using namespace std;

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

void extractMesh(Handle_Geom_BSplineSurface, std::vector<float>&, int);

void writeMesh(std::vector<float>&, int, int, std::string);

void writeHexahedralMesh(std::vector<float>& , int , int , int , std::string);

#endif