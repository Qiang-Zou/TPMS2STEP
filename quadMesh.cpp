#include "quadMesh.h"

void extractMesh(Handle_Geom_BSplineSurface patch, std::vector<float>& points, int numberOfMeshPoints) {
    for (int i=0;i<numberOfMeshPoints;i++) {
        for (int j=0;j<numberOfMeshPoints;j++) {
            float u = j*1.0/(numberOfMeshPoints-1);
            float v = i*1.0/(numberOfMeshPoints-1);
            gp_Pnt point;
            patch->D0(u,v,point);
            points.push_back(point.X());
            points.push_back(point.Y());
            points.push_back(point.Z());
        }
    }
    return;
}

void writeMesh(std::vector<float>& points, int numberOfFaces, int numberOfMeshPoints, std::string fileName) {
    MyMesh mesh;
    std::vector<MyMesh::VertexHandle> tmp_face_vhandles;
    std::vector<MyMesh::VertexHandle> vertices;
    for (int m=0;m<numberOfFaces;m++) {
        vertices.clear();
        for (int i=0;i<numberOfMeshPoints;i++) {
            for (int j=0;j<numberOfMeshPoints;j++) {
                float x = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3];
                float y = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3+1];
                float z = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3+2];
                MyMesh::VertexHandle vhandle = mesh.add_vertex(MyMesh::Point(x,y,z));
                vertices.push_back(vhandle);
            }
        }
        for (int i=0;i<numberOfMeshPoints-1;i++) {
            for (int j=0;j<numberOfMeshPoints-1;j++) {
                tmp_face_vhandles.clear();
                tmp_face_vhandles.push_back(vertices[i*numberOfMeshPoints+j]);
                tmp_face_vhandles.push_back(vertices[i*numberOfMeshPoints+j+1]);
                tmp_face_vhandles.push_back(vertices[(i+1)*numberOfMeshPoints+(j+1)]);
                tmp_face_vhandles.push_back(vertices[(i+1)*numberOfMeshPoints+j]);
                mesh.add_face(tmp_face_vhandles);
            }
        }
    }
    try {
        if (!OpenMesh::IO::write_mesh(mesh, fileName)) {
            std::cerr << "TPMS2STEP > Write mesh error!\n";
            return;
        }
    }
    catch (std::exception& x) {
        std::cerr << x.what() << std::endl;
        return;
    }
    cout << "TPMS2STEP > Write mesh ok!" << endl;
}

void writeHexahedralMesh(std::vector<float>& points, int numberOfFaces, int numberOfFacesPerCell, int numberOfMeshPoints, std::string fileName) {
    MyMesh mesh_1;
    MyMesh mesh_2;
    std::vector<MyMesh::VertexHandle> tmp_face_vhandles;
    std::vector<MyMesh::VertexHandle> vertices;
    for (int r=0;r<numberOfFaces/numberOfFacesPerCell;r+=2) {
        for (int s=0;s<numberOfFacesPerCell;s++) {
            int m = r*numberOfFacesPerCell+s;
            vertices.clear();
            for (int i=0;i<numberOfMeshPoints;i++) {
                for (int j=0;j<numberOfMeshPoints;j++) {
                    float x = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3];
                    float y = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3+1];
                    float z = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3+2];
                    MyMesh::VertexHandle vhandle = mesh_1.add_vertex(MyMesh::Point(x,y,z));
                    vertices.push_back(vhandle);
                }
            }
            for (int i=0;i<numberOfMeshPoints-1;i++) {
                for (int j=0;j<numberOfMeshPoints-1;j++) {
                    // cout << "i: " << i << "  j: " << j << endl;
                    tmp_face_vhandles.clear();
                    tmp_face_vhandles.push_back(vertices[i*numberOfMeshPoints+j]);
                    tmp_face_vhandles.push_back(vertices[i*numberOfMeshPoints+j+1]);
                    tmp_face_vhandles.push_back(vertices[(i+1)*numberOfMeshPoints+(j+1)]);
                    tmp_face_vhandles.push_back(vertices[(i+1)*numberOfMeshPoints+j]);
                    mesh_1.add_face(tmp_face_vhandles);
                }
            }
        }
    }
    for (int r=1;r<numberOfFaces/numberOfFacesPerCell;r+=2) {
        for (int s=0;s<numberOfFacesPerCell;s++) {
            int m = r*numberOfFacesPerCell+s;
            vertices.clear();
            for (int i=0;i<numberOfMeshPoints;i++) {
                for (int j=0;j<numberOfMeshPoints;j++) {
                    float x = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3];
                    float y = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3+1];
                    float z = points[(m*numberOfMeshPoints*numberOfMeshPoints+(i*numberOfMeshPoints+j))*3+2];
                    MyMesh::VertexHandle vhandle = mesh_2.add_vertex(MyMesh::Point(x,y,z));
                    vertices.push_back(vhandle);
                }
            }
            for (int i=0;i<numberOfMeshPoints-1;i++) {
                for (int j=0;j<numberOfMeshPoints-1;j++) {
                    // cout << "i: " << i << "  j: " << j << endl;
                    tmp_face_vhandles.clear();
                    tmp_face_vhandles.push_back(vertices[i*numberOfMeshPoints+j]);
                    tmp_face_vhandles.push_back(vertices[i*numberOfMeshPoints+j+1]);
                    tmp_face_vhandles.push_back(vertices[(i+1)*numberOfMeshPoints+(j+1)]);
                    tmp_face_vhandles.push_back(vertices[(i+1)*numberOfMeshPoints+j]);
                    mesh_2.add_face(tmp_face_vhandles);
                }
            }
        }
    }
    // get the vertices of each quad mesh
    std::vector<MyMesh::Point> verticesHex_1;
    std::vector<MyMesh::Point> verticesHex_2;
    for (MyMesh::FaceIter f_it=mesh_1.faces_begin(); f_it!=mesh_1.faces_end(); ++f_it) {
        for (MyMesh::FaceVertexIter fv_it=mesh_1.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
            MyMesh::Point tempPoint = mesh_1.point(*fv_it);
            verticesHex_1.push_back(tempPoint);
        }
    }
    for (MyMesh::FaceIter f_it=mesh_2.faces_begin(); f_it!=mesh_2.faces_end(); ++f_it) {
        for (MyMesh::FaceVertexIter fv_it=mesh_2.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
            MyMesh::Point tempPoint = mesh_2.point(*fv_it);
            verticesHex_2.push_back(tempPoint);
        }
    }
    // write to file
    ofstream ofs;
    ofs.open(fileName, ios::out);
    ofs << "OFF\n" << to_string(numberOfFaces*(numberOfMeshPoints-1)*(numberOfMeshPoints-1)*4) << " 0 0" << endl;
    for (int i=0;i<(numberOfFaces/2);i++) {
        for (int j=0;j<(numberOfMeshPoints-1)*(numberOfMeshPoints-1);j++) {
            for (int k=0;k<4;k++) {
                ofs << verticesHex_1[i*((numberOfMeshPoints-1)*(numberOfMeshPoints-1)*4)+j*4+k][0] << " "
                << verticesHex_1[i*((numberOfMeshPoints-1)*(numberOfMeshPoints-1)*4)+j*4+k][1] << " "
                << verticesHex_1[i*((numberOfMeshPoints-1)*(numberOfMeshPoints-1)*4)+j*4+k][2] << endl;
            }
            for (int k=0;k<4;k++) {
                ofs << verticesHex_2[i*((numberOfMeshPoints-1)*(numberOfMeshPoints-1)*4)+j*4+k][0] << " "
                << verticesHex_2[i*((numberOfMeshPoints-1)*(numberOfMeshPoints-1)*4)+j*4+k][1] << " "
                << verticesHex_2[i*((numberOfMeshPoints-1)*(numberOfMeshPoints-1)*4)+j*4+k][2] << endl;
            }
        }
    }
    ofs.close();
    cout << "TPMS2STEP > Write hexahedral mesh ok!" << endl;
}
