#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

#include "Polyhedron.hpp"

using namespace std;

using namespace Eigen;
namespace PolyhedronLibrary
{

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
        for (auto& v : verts) v = v.normalized();
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
        for (auto& v : verts) v = v.normalized();
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

        for (auto& v : raw) verts.push_back(v.normalized());

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
Polyhedron buildPlatonicSolid(unsigned int& q) 
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

        const Vector3d v = verts[i];
        //ector<double> coords = {v(0), v(1), v(2)};
        //P.IdCell0Ds[i] = coords; 
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
        for (unsigned int i = 0; i < f.size(); ++i)
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
// Funzione per la triangolazione di classe 1
Polyhedron triangulateClass1(PolyhedronLibrary::Polyhedron& P, unsigned int& t_value)
{
    Polyhedron P_triangolato;
    
    // Copia dei dati iniziali
    unsigned int nextPointId = 0; // nextPointId è inizializzata con il numero inizale di vertici
    unsigned int nextEdgeId = 0;
    unsigned int nextFaceId = 0;

    map<pair<unsigned int, unsigned int>, unsigned int> edgeMap; // mappa che tiene traccia dei lati già esistenti evitando duplicati
    map<tuple<double, double, double>, unsigned int> pointMap;

    // LATI
    vector<vector<unsigned int>> originalFaces = P.Cell2DsVertices;
    // per ogni faccia salvo gli id vertici A B C
    for (size_t fid = 0; fid < originalFaces.size(); ++fid) 
    {
        const auto& face = P.Cell2DsVertices[fid];

        unsigned int A = face[0]; 
        unsigned int B = face[1];
        unsigned int C = face[2];

        // salvo le tre coordinate di ogni vertice
        Vector3d vA = P.Cell0DsCoordinates.col(A);
        Vector3d vB = P.Cell0DsCoordinates.col(B);
        Vector3d vC = P.Cell0DsCoordinates.col(C);

        vector<vector<unsigned int>> rows;

        // Costruisci i punti interni con interpolazione baricentrica
        for (unsigned int i = 0; i < t_value+1; ++i) 
        { // ogni iterazione rappresenta una riga di punti nel triangolo
            vector<unsigned int> row; // contiene gli id dei punti generati in quella riga
            // ogni iterazione rappresenta un punto sulla riga
            for (unsigned int j = 0; j <= i; ++j) 
            { 
                // trovo i parametri per dividere il lato i t_value parti
                double alpha = double(t_value - i) / t_value; 
                double beta  = double(i - j) / t_value;
                double gamma = double(j) / t_value;
                Vector3d point = alpha * vA + beta * vB + gamma * vC; // coordinate del nuovo punto
                //point = point.normalized(); // proiezione dei vertici sulla sfera con la normalizzazione

                // Verifica duplicati 
                auto key = make_tuple(point(0), point(1), point(2));
                auto it = pointMap.find(key);
                if (it != pointMap.end()) 
                {
                    // Punto già esistente
                    row.push_back(it->second);
                } 
                else 
                {
                    // Nuovo punto
                    P_triangolato.Cell0DsCoordinates.conservativeResize(3, nextPointId + 1);
                    P_triangolato.Cell0DsCoordinates.col(nextPointId) = point;
                    P_triangolato.Cell0DsId.push_back(nextPointId);
                    //P_triangolato.IdCell0Ds[nextPointId] = {point(0), point(1), point(2)};
                    row.push_back(nextPointId);
                    pointMap[key] = nextPointId;
                    ++nextPointId;
                }  
            }
            rows.push_back(row);
        }

        // Ora crea triangoli
        for (unsigned int i = 0; i < t_value; ++i) 
        { //ogni iterazione rappresenta una riga di triangoli che verranno creati
            for (unsigned int j = 0; j < i; ++j) 
            { //ogni iterazione rappresenta un triangolo che verrà creato nella riga corrente

                // triangolazione di tutta la faccia tranne il triangolo in cima, a partire dalla linea con due punti
                vector<unsigned int> tri1 = { rows[i][j], rows[i + 1][j], rows[i + 1][j + 1] };
                vector<unsigned int> tri2 = { rows[i][j], rows[i][j + 1], rows[i + 1][j + 1] };

                for (auto& tri : {tri1, tri2}) 
                {
                    P_triangolato.Cell2DsVertices.push_back(tri);
                    P_triangolato.Cell2DsId.push_back(nextFaceId++);

                    vector<unsigned int> edgeIds; //lista degli id dei lati del triangolo corrente
                    for (int k = 0; k < 3; ++k) 

                    
                    { // itera su igni vertice del triangolo corrente
                        auto v1 = tri[k]; // vertice corrente
                        auto v2 = tri[(k + 1) % 3]; // vertice successivo + condizione per la chiusura del tiangolo
                        auto key = minmax(v1, v2); // ordina la coppia per evitare duplicati
                        if (edgeMap.find(key) == edgeMap.end()) 
                        {
                            edgeMap[key] = nextEdgeId++;
                            P_triangolato.Cell1DsId.push_back(edgeMap[key]); // aggiunge l'id del nuovo lato a Cell1DsId
                        }
                        edgeIds.push_back(edgeMap[key]);
                    }
                    P_triangolato.Cell2DsEdges.push_back(edgeIds);
                }
            }

            // triangolazione del triangolo in cima
            vector<unsigned int> tri = { rows[i][i], rows[i + 1][i], rows[i + 1][i + 1] };
            P_triangolato.Cell2DsVertices.push_back(tri);
            P_triangolato.Cell2DsId.push_back(nextFaceId++);

            vector<unsigned int> edgeIds;
            for (int k = 0; k < 3; ++k) {
                auto v1 = tri[k];
                auto v2 = tri[(k + 1) % 3];
                auto key = minmax(v1, v2);
                if (edgeMap.find(key) == edgeMap.end()) 
                {
                    edgeMap[key] = nextEdgeId++;
                    P_triangolato.Cell1DsId.push_back(edgeMap[key]);
                }
                edgeIds.push_back(edgeMap[key]);
            }
            P_triangolato.Cell2DsEdges.push_back(edgeIds);
        }
    }

    // Aggiorniamo gli attributi del poliedro
    P_triangolato.NumCell0Ds = nextPointId;
    P_triangolato.NumCell1Ds = edgeMap.size();
    P_triangolato.NumCell2Ds = nextFaceId;

    P_triangolato.Cell1DsExtrema = MatrixXi(2, P_triangolato.NumCell1Ds);
    // viene riempita la matrice Cell1DsExtrema con i valori aggiornati di edgeMap
    for (const auto& [key, eid] : edgeMap) 
    {
        P_triangolato.Cell1DsExtrema(0, eid) = key.first;
        P_triangolato.Cell1DsExtrema(1, eid) = key.second;
    }

    // Aggiornamento Cell3D
    P_triangolato.Cell3DsId = {0};
    P_triangolato.Cell3DsVertices = {P_triangolato.Cell0DsId};
    P_triangolato.Cell3DsEdges = {P_triangolato.Cell1DsId};
    P_triangolato.Cell3DsFaces = {P_triangolato.Cell2DsId};
    P_triangolato.NumCell3Ds = 1;

    return P_triangolato;
}

/*
// ***************************************************************************
// Funzione per la triangolazione di classe 2
Polyhedron triangulateClass2(PolyhedronLibrary::Polyhedron& P, unsigned int& t_value)
{
    using namespace Eigen;
    Polyhedron P_triangolato1 = triangulateClass1(P, t_value); //triangolazione di classe 1
    Polyhedron P_triangolato2; //creo nuova mesh per salvare dati della nuova traingolazione

    unsigned int nextPointId = 0;
    unsigned int nextEdgeId = 0;
    unsigned int nextFaceId = 0;

    map<tuple<double, double, double>, unsigned int> pointMap;
    map<pair<unsigned int, unsigned int>, unsigned int> edgeMap;
    //unordered_set<unsigned int> midpoints_ids;
    set<unsigned int> midpoints_ids;

    for (const auto& face : P_triangolato1.Cell2DsVertices)  //itera su ogni faccia della triangolazione di classe 1 
    {
        unsigned int A = face[0];
        unsigned int B = face[1];
        unsigned int C = face[2];

        Vector3d vA = P_triangolato1.Cell0DsCoordinates.col(A);
        Vector3d vB = P_triangolato1.Cell0DsCoordinates.col(B);
        Vector3d vC = P_triangolato1.Cell0DsCoordinates.col(C);

        Vector3d M_ab = 0.5 * (vA + vB); //punto medio AB
        M_ab = M_ab.normalized();
        Vector3d M_bc = 0.5 * (vB + vC);
        M_bc = M_bc.normalized();
        Vector3d M_ca = 0.5 * (vC + vA);
        M_ca = M_ca.normalized();

        Vector3d centroid = (vA + vB + vC) / 3; //baricentro
        centroid = centroid.normalized();

        vector<vector<unsigned int>> rows;

                // Verifica duplicati 
                auto key = make_tuple(point(0), point(1), point(2));
                auto it = pointMap.find(key);
                if (it != pointMap.end()) 
                {
                    // Punto già esistente
                    row.push_back(it->second);
                } 
                else 
                {
                    // Nuovo punto
                    P_triangolato.Cell0DsCoordinates.conservativeResize(3, nextPointId + 1);
                    P_triangolato.Cell0DsCoordinates.col(nextPointId) = point;
                    P_triangolato.Cell0DsId.push_back(nextPointId);
                    P_triangolato.IdCell0Ds[nextPointId] = {point(0), point(1), point(2)};
                    row.push_back(nextPointId);
                    pointMap[key] = nextPointId;
                    ++nextPointId;
                }  
            }
            rows.push_back(row);
        }

        for (unsigned int i = 0; i < t_value; ++i)
        {
            for (unsigned int j = 0; j < i; ++j)
            {
                std::vector<unsigned int> tri1 = { rows[i][j], rows[i + 1][j], rows[i + 1][j + 1] };
                std::vector<unsigned int> tri2 = { rows[i][j], rows[i][j + 1], rows[i + 1][j + 1] };

                std::vector<std::vector<unsigned int>> triangles = { tri1, tri2 };

                for (const auto& tri : triangles)
                {
                    P_triangolato.Cell2DsVertices.push_back(tri);
                    P_triangolato.Cell2DsId.push_back(nextFaceId++);

                    std::vector<unsigned int> edgeIds;
                    for (int k = 0; k < 3; ++k)
                    {
                        unsigned int v1 = tri[k];
                        unsigned int v2 = tri[(k + 1) % 3];
                        auto key = std::minmax(v1, v2);
                        if (edgeMap.find(key) == edgeMap.end())
                        {
                            edgeMap[key] = nextEdgeId++;
                            P_triangolato.Cell1DsId.push_back(edgeMap[key]);
                        }
                        edgeIds.push_back(edgeMap[key]);
                    }
                    P_triangolato.Cell2DsEdges.push_back(edgeIds);
                }
            }

            // Triangolo sul bordo diagonale
            vector<unsigned int> tri = { rows[i][i], rows[i + 1][i], rows[i + 1][i + 1] };
            P_triangolato.Cell2DsVertices.push_back(tri);
            P_triangolato.Cell2DsId.push_back(nextFaceId++);

            vector<unsigned int> edgeIds;
            for (int k = 0; k < 3; ++k)
            {
                unsigned int v1 = tri[k];
                unsigned int v2 = tri[(k + 1) % 3];
                auto key = std::minmax(v1, v2);
                if (edgeMap.find(key) == edgeMap.end())
                {
                    edgeMap[key] = nextEdgeId++;
                    P_triangolato.Cell1DsId.push_back(edgeMap[key]);
                }
                edgeIds.push_back(edgeMap[key]);
            }
            P_triangolato.Cell2DsEdges.push_back(edgeIds);
        }
    }

    P_triangolato.NumCell0Ds = nextPointId;
    P_triangolato.NumCell1Ds = edgeMap.size();
    P_triangolato.NumCell2Ds = nextFaceId;

    P_triangolato.Cell1DsExtrema = MatrixXi(2, P_triangolato.NumCell1Ds);
    for (const auto& [key, eid] : edgeMap)
    {
        P_triangolato.Cell1DsExtrema(0, eid) = key.first;
        P_triangolato.Cell1DsExtrema(1, eid) = key.second;
    }

    P_triangolato.Cell3DsId = { 0 };
    P_triangolato.Cell3DsVertices = { P_triangolato.Cell0DsId };
    P_triangolato.Cell3DsEdges = { P_triangolato.Cell1DsId };
    P_triangolato.Cell3DsFaces = { P_triangolato.Cell2DsId };
    P_triangolato.NumCell3Ds = 1;

    for (const auto& [id, coords] : P_triangolato.IdCell0Ds) {
        std::cout << "ID: " << id << " -> [";
        for (size_t i = 0; i < coords.size(); ++i) {
            std::cout << coords[i];
            if (i < coords.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }

    return P_triangolato;
}
*/

// ***************************************************************************
// Funzione per la dualizzazione
Polyhedron Dualize(const Polyhedron& P_original) 
{
    Polyhedron P_duale;
    // 1. Ogni faccia originale diventa un vertice del duale
    unsigned int num_faces = P_original.Cell2DsVertices.size();
    P_duale.NumCell0Ds = num_faces;
    P_duale.Cell0DsId.resize(num_faces);
    P_duale.Cell0DsCoordinates.resize(3, num_faces);

    for (unsigned int i = 0; i < num_faces; ++i) {
        P_duale.Cell0DsId[i] = i;

        const vector<unsigned int> face = P_original.Cell2DsVertices[i];
        Vector3d centroid(0, 0, 0);
        for (unsigned int vid : face) 
        {
            centroid += P_original.Cell0DsCoordinates.col(vid);
        }
        centroid = centroid / face.size();
        centroid.normalize(); // Proiezione sulla sfera
        P_duale.Cell0DsCoordinates.col(i) = centroid;
    }

    // 2. Costruisci mappa edge -> facce (per trovare adiacenze tra facce)
    std::map<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>> edge_to_faces;
    for (unsigned int i = 0; i < num_faces; ++i) {
        const std::vector<unsigned int> face = P_original.Cell2DsVertices[i];
        unsigned int n = face.size();
        for (unsigned int j = 0; j < n; ++j) {
            unsigned int v1 = face[j];
            unsigned int v2 = face[(j + 1) % n];
            auto edge = std::minmax(v1, v2);
            edge_to_faces[edge].push_back(i);
        }
    }

    // 3. Ogni coppia di facce adiacenti genera un lato nel duale
    unsigned int edge_id = 0;
    P_duale.Cell1DsExtrema.resize(2, edge_to_faces.size());
    for (const auto& [edge, faces] : edge_to_faces) 
    {
        if (faces.size() == 2) 
        {
            unsigned int f1 = faces[0];
            unsigned int f2 = faces[1];
            P_duale.Cell1DsId.push_back(edge_id);
            P_duale.Cell1DsExtrema.col(edge_id) << f1, f2;
        }
    }

    P_duale.NumCell1Ds = edge_id;
    P_duale.Cell1DsExtrema.conservativeResize(2, edge_id);

    // 4. Ogni vertice originale diventa una faccia nel duale
    std::map<unsigned int, std::vector<unsigned int>> vertex_to_faces;
    for (unsigned int i = 0; i < num_faces; ++i) {
        for (unsigned int vid : P_original.Cell2DsVertices[i]) {
            vertex_to_faces[vid].push_back(i);
        }
    }

    unsigned int face_id = 0;
    for (const auto& [vid, faces] : vertex_to_faces) {
        std::vector<unsigned int> face = faces;

        P_duale.Cell2DsId.push_back(face_id);
        P_duale.Cell2DsVertices.push_back(face);
        ++face_id;
    }
    P_duale.NumCell2Ds = face_id;

    // 5. Assegna l’unica cella 3D del duale
    P_duale.NumCell3Ds = 1;
    P_duale.Cell3DsId = {0};
    P_duale.Cell3DsFaces.resize(1);
    for (unsigned int i = 0; i < P_duale.NumCell2Ds; ++i) {
        P_duale.Cell3DsFaces[0].push_back(i);
    }

    // Aggiornamento Cell3D
    P_duale.Cell3DsId = {0};
    P_duale.Cell3DsVertices = {P_duale.Cell0DsId};
    P_duale.Cell3DsEdges = {P_duale.Cell1DsId};
    P_duale.Cell3DsFaces = {P_duale.Cell2DsId};
    P_duale.NumCell3Ds = 1;

    // Riempimento Cell2DsEdges
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> edge_lookup;
    for (unsigned int eid = 0; eid < P_duale.NumCell1Ds; ++eid) {
        unsigned int v1 = P_duale.Cell1DsExtrema(0, eid);
        unsigned int v2 = P_duale.Cell1DsExtrema(1, eid);
        auto key = std::minmax(v1, v2);
        edge_lookup[key] = eid;
    }

    P_duale.Cell2DsEdges.resize(P_duale.NumCell2Ds);
    for (unsigned int fid = 0; fid < P_duale.NumCell2Ds; ++fid) 
    {
        const std::vector<unsigned int> verts = P_duale.Cell2DsVertices[fid];
        std::vector<unsigned int> edge_ids;
        unsigned int n = verts.size();

        for (unsigned int i = 0; i < n; ++i) 
        {
            unsigned int v1 = verts[i];
            unsigned int v2 = verts[(i + 1) % n];
            auto key = std::minmax(v1, v2);

            auto it = edge_lookup.find(key);
            if (it != edge_lookup.end()) {
                edge_ids.push_back(it->second);
            } 
            else 
            {   
                edge_lookup[key] = edge_id++;
            }
            edge_ids.push_back(edge_lookup[key]);
        }
        P_duale.Cell2DsEdges[fid] = edge_ids;
    }
}

}