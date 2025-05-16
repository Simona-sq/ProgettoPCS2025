#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

#include "Polyhedron.hpp"

namespace PolyhedronLibrary
{
// ***************************************************************************
// Funzione di normalizzazione 
Vector3d normalize(Vector3d v) 
{
    return v.normalized(); // vettore normalizzato
}

// ***************************************************************************
// Definisco i vertici e le facce del solido platonico di partenza
pair<vector<Vector3d>, vector<vector<unsigned int>>> getSolidData(unsigned int& q) 
{
    vector<Vector3d> verts;
    vector<vector<unsigned int>> faces;

    if (q == 3) {
        // Tetraedro regolare
        verts = {
            {0.0, 0.0, sqrt(3)},
            {sqrt(8.0/3.0), 0.0, -1.0/sqrt(3)},
            {-sqrt(2.0/3.0), sqrt(2.0), -1.0/sqrt(3)},
            {-sqrt(2.0/3.0), -sqrt(2.0), -1.0/sqrt(3)}
        };
        for (auto& v : verts) v = normalize(v);
        faces = {
            {0, 1, 2},
            {0, 2, 3},
            {0, 3, 1},
            {1, 3, 2}
        };
    } else if (q == 4) {
        // Ottaedro
        verts = {
            {1, 0, 0}, {-1, 0, 0},
            {0, 1, 0}, {0, -1, 0},
            {0, 0, 1}, {0, 0, -1}
        };
        for (auto& v : verts) v = normalize(v);
        faces = {
            {0, 4, 2}, {2, 4, 1}, {1, 4, 3}, {3, 4, 0},
            {0, 2, 5}, {2, 1, 5}, {1, 3, 5}, {3, 0, 5}
        };
    } else if (q == 5) {
        // Icosaedro
        vector<Vector3d> raw = {
            { 0.0, 0.0, 1.175571}, { 1.051462, 0.0, 0.5257311}, {0.3249197, 1.0, 0.5257311},
            { -0.8506508, 0.618034, 0.5257311}, {-0.8506508, -0.618034, 0.5257311}, { 0.3249197, -1.0, 0.5257311}, 
            { 0.8506508, 0.618034, -0.5257311}, { 0.8506508, -0.618034, -0.5257311}, {-0.3249197, 1.0, -0.5257311}, 
            {-1.051462, 0.0, -0.5257311}, {-0.3249197, -1.0, -0.5257311}, {0.0, 0.0, -1.175571}
        };

        for (auto& v : raw) verts.push_back(normalize(v));

        faces = {
            {0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 5, 1},
            {1, 5, 7}, {1, 7, 6}, {1, 6, 2}, {2, 6, 8}, {2, 8, 3},
            {3, 8, 9}, {3, 9, 4}, {4, 9, 10}, {4, 10, 5}, {5, 10, 7},
            {6, 7, 11}, {6, 11, 8}, {7, 10, 11}, {8, 11, 9}, {9, 11, 10}
        };
    }

    return {verts, faces};
}

// ***************************************************************************
// Riempio la struttura Polyhedron
Polyhedron buildPlatonicSolid(unsigned int& p, unsigned int& q, unsigned int& b, unsigned int& c) 
{
    Polyhedron P;

    auto [verts, faces] = getSolidData(q);

    // Riempimento Cell0Ds
    P.NumCell0Ds = verts.size();
    P.Cell0DsCoordinates = MatrixXd(3, P.NumCell0Ds);
    for (unsigned int i = 0; i < verts.size(); ++i) 
    {
        P.Cell0DsId.push_back(i);
        P.Cell0DsCoordinates.col(i) = verts[i];
        P.IdCell0Ds , verts , verts;
    }
    

    // Creazione lati univoci con orientamento coerente
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> edgeMap;
    unsigned int edgeCounter = 0;

    for (unsigned int fid = 0; fid < faces.size(); fid++) 
    {
        const auto& f = faces[fid];
        P.Cell2DsId.push_back(fid);
        P.Cell2DsVertices.push_back(f);

        std::vector<unsigned int> edgeIds;
        for (int i = 0; i < f.size(); ++i)
        {
            unsigned int v1 = f[i];
            unsigned int v2 = f[(i + 1) % f.size()];
            std::pair<unsigned int, unsigned int> key = {v1, v2};
            std::pair<unsigned int, unsigned int> revKey = {v2, v1};

            if (edgeMap.find(key) == edgeMap.end() && edgeMap.find(revKey) == edgeMap.end()) 
            {
                edgeMap[key] = edgeCounter;
                P.Cell1DsId.push_back(edgeCounter);
                edgeCounter++;
            }

            if (edgeMap.find(key) != edgeMap.end()) 
            {
                edgeIds.push_back(edgeMap[key]);
            } 
            else 
            {
                edgeIds.push_back(edgeMap[revKey]);
            }
        }
        P.Cell2DsEdges.push_back(edgeIds);
    }

    // Riempimento Cell1Ds coerente con chiave usata
    P.NumCell1Ds = edgeMap.size();
    P.Cell1DsExtrema = MatrixXi(2, P.NumCell1Ds);
    for (const auto& [key, eid] : edgeMap) 
    {
        P.Cell1DsExtrema(0, eid) = key.first;
        P.Cell1DsExtrema(1, eid) = key.second;
    }

    // Cell3Ds (1 solo poliedro)
    P.NumCell2Ds = faces.size();
    P.NumCell3Ds = 1;
    P.Cell3DsId.push_back(0);
    P.Cell3DsVertices.push_back(P.Cell0DsId);
    P.Cell3DsEdges.push_back(P.Cell1DsId);
    P.Cell3DsFaces.push_back(P.Cell2DsId);

    return P;
}


// ***************************************************************************
// Funzione per la triangolazione
void triangulateClass1(PolyhedronLibrary::Polyhedron& P, unsigned int& t_value)
{
    using namespace Eigen;

    // Copia dei dati iniziali
    unsigned int nextPointId = P.NumCell0Ds;
    unsigned int nextEdgeId = P.NumCell1Ds;
    unsigned int nextFaceId = P.NumCell2Ds;

    std::map<std::pair<unsigned int, unsigned int>, unsigned int> edgeMap;

    // Prepara i vecchi spigoli nel map per evitare duplicati
    for (unsigned int eid = 0; eid < P.NumCell1Ds; ++eid) {
        auto v1 = P.Cell1DsExtrema(0, eid);
        auto v2 = P.Cell1DsExtrema(1, eid);
        auto key = std::minmax(v1, v2);
        edgeMap[key] = eid;
    }

    for (size_t fid = 0; fid < P.Cell2DsVertices.size(); ++fid) {
        const auto& face = P.Cell2DsVertices[fid];
        if (face.size() != 3) {
            std::cerr << "Solo triangoli supportati per la classe I\n";
            continue;
        }

        unsigned int A = face[0];
        unsigned int B = face[1];
        unsigned int C = face[2];

        Vector3d vA = P.Cell0DsCoordinates.col(A);
        Vector3d vB = P.Cell0DsCoordinates.col(B);
        Vector3d vC = P.Cell0DsCoordinates.col(C);

        std::vector<std::vector<unsigned int>> rows;

        // Costruisci i punti interni con interpolazione baricentrica
        for (unsigned int i = 0; i <= t_value; ++i) {
            std::vector<unsigned int> row;
            for (unsigned int j = 0; j <= i; ++j) {
                double alpha = double(t_value - i) / t_value;
                double beta  = double(i - j) / t_value;
                double gamma = double(j) / t_value;
                Vector3d point = alpha * vA + beta * vB + gamma * vC;

                // Verifica duplicati (facoltativo)
                P.Cell0DsCoordinates.conservativeResize(3, nextPointId + 1);
                P.Cell0DsCoordinates.col(nextPointId) = point;
                P.Cell0DsId.push_back(nextPointId);
                P.IdCell0Ds[nextPointId] = {point(0), point(1), point(2)};
                row.push_back(nextPointId);
                ++nextPointId;
            }
            rows.push_back(row);
        }

        // Ora crea triangoli
        for (unsigned int i = 0; i < t_value; ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                // triangolo in basso
                std::vector<unsigned int> tri1 = { rows[i][j], rows[i + 1][j], rows[i + 1][j + 1] };
                std::vector<unsigned int> tri2 = { rows[i][j], rows[i][j + 1], rows[i + 1][j + 1] };

                for (auto& tri : {tri1, tri2}) {
                    P.Cell2DsVertices.push_back(tri);
                    P.Cell2DsId.push_back(nextFaceId++);

                    std::vector<unsigned int> edgeIds;
                    for (int k = 0; k < 3; ++k) {
                        auto v1 = tri[k];
                        auto v2 = tri[(k + 1) % 3];
                        auto key = std::minmax(v1, v2);
                        if (edgeMap.find(key) == edgeMap.end()) {
                            edgeMap[key] = nextEdgeId++;
                            P.Cell1DsId.push_back(edgeMap[key]);
                        }
                        edgeIds.push_back(edgeMap[key]);
                    }
                    P.Cell2DsEdges.push_back(edgeIds);
                }
            }

            // Triangolo in alto (ultima colonna)
            std::vector<unsigned int> tri = { rows[i][i], rows[i + 1][i], rows[i + 1][i + 1] };
            P.Cell2DsVertices.push_back(tri);
            P.Cell2DsId.push_back(nextFaceId++);

            std::vector<unsigned int> edgeIds;
            for (int k = 0; k < 3; ++k) {
                auto v1 = tri[k];
                auto v2 = tri[(k + 1) % 3];
                auto key = std::minmax(v1, v2);
                if (edgeMap.find(key) == edgeMap.end()) {
                    edgeMap[key] = nextEdgeId++;
                    P.Cell1DsId.push_back(edgeMap[key]);
                }
                edgeIds.push_back(edgeMap[key]);
            }
            P.Cell2DsEdges.push_back(edgeIds);
        }
    }

    // Riempimento Cell1DsExtrema finale
    P.NumCell0Ds = nextPointId;
    P.NumCell1Ds = edgeMap.size();
    P.NumCell2Ds = nextFaceId;

    P.Cell1DsExtrema = MatrixXi(2, P.NumCell1Ds);
    for (const auto& [key, eid] : edgeMap) {
        P.Cell1DsExtrema(0, eid) = key.first;
        P.Cell1DsExtrema(1, eid) = key.second;
    }

    // Aggiornamento Cell3D
    P.Cell3DsId = {0};
    P.Cell3DsVertices = {P.Cell0DsId};
    P.Cell3DsEdges = {P.Cell1DsId};
    P.Cell3DsFaces = {P.Cell2DsId};
    P.NumCell3Ds = 1;

    //return P;
}
}
