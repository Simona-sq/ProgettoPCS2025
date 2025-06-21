#include "Utils.hpp"
#include "UCDUtilities.hpp"
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <set>

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
    for (unsigned int i = 0; i < verts.size(); i++) 
    {
        P.Cell0DsId.push_back(i);
        P.Cell0DsCoordinates.col(i) = verts[i];
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
        for (unsigned int i = 0; i < f.size(); i++)
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
    P.Cell3DsVertices.push_back(P.Cell0DsId);
    P.Cell3DsEdges.push_back(P.Cell1DsId);
    P.Cell3DsFaces.push_back(P.Cell2DsId);

    return P;
}


// ***************************************************************************
// Funzione per la triangolazione di classe 1
Polyhedron triangulateClass1(const Polyhedron& P, const unsigned int& t_value)
{
    Polyhedron P_triangolato;
    
    // Copia dei dati iniziali
    unsigned int nextPointId = 0;
    unsigned int nextEdgeId = 0;
    unsigned int nextFaceId = 0;

    map<pair<unsigned int, unsigned int>, unsigned int> edgeMap; // mappa che tiene traccia dei lati già esistenti evitando duplicati
    map<tuple<double, double, double>, unsigned int> pointMap;

    // LATI
    vector<vector<unsigned int>> originalFaces = P.Cell2DsVertices;
    
    for (size_t fid = 0; fid < originalFaces.size(); fid++) 
    {
        const auto& face = P.Cell2DsVertices[fid];

        // per ogni faccia salvo gli id vertici A B C
        unsigned int A = face[0]; 
        unsigned int B = face[1];
        unsigned int C = face[2];

        // salvo le tre coordinate di ogni vertice
        Vector3d vA = P.Cell0DsCoordinates.col(A);
        Vector3d vB = P.Cell0DsCoordinates.col(B);
        Vector3d vC = P.Cell0DsCoordinates.col(C);

        vector<vector<unsigned int>> rows;

        // Costruisci i punti interni con interpolazione baricentrica
        for (unsigned int i = 0; i < t_value+1; i++) 
        { // ogni iterazione rappresenta una riga di punti nel triangolo
            vector<unsigned int> row; // contiene gli id dei punti generati in quella riga
            // ogni iterazione rappresenta un punto sulla riga
            for (unsigned int j = 0; j <= i; j++) 
            { 
                // trovo i parametri per dividere il lato i t_value parti
                double alpha = double(t_value - i) / t_value; 
                double beta  = double(i - j) / t_value;
                double gamma = double(j) / t_value;
                Vector3d point = alpha * vA + beta * vB + gamma * vC; // coordinate del nuovo punto
                
                // Verifica duplicati 
                auto key = make_tuple(point(0), point(1), point(2));
                auto it = pointMap.find(key);
                if (it != pointMap.end()) //il punto esiste
                {
                    row.push_back(it->second);
                } 
                else //il punto è nuovo
                {
                    // Nuovo punto
                    P_triangolato.Cell0DsCoordinates.conservativeResize(3, nextPointId + 1);
                    P_triangolato.Cell0DsCoordinates.col(nextPointId) = point;
                    P_triangolato.Cell0DsId.push_back(nextPointId);
                    row.push_back(nextPointId);
                    pointMap[key] = nextPointId;
                    ++nextPointId;
                }  
            }
            rows.push_back(row);
        }

        // Ora crea triangoli
        for (unsigned int i = 0; i < t_value; i++) 
        { //ogni iterazione rappresenta una riga di triangoli che verranno creati
            for (unsigned int j = 0; j < i; j++) 
            { //ogni iterazione rappresenta un triangolo che verrà creato nella riga corrente

                // triangolazione di tutta la faccia tranne il triangolo in cima, a partire dalla linea con due punti
                vector<unsigned int> tri1 = { rows[i][j], rows[i + 1][j], rows[i + 1][j + 1] };
                vector<unsigned int> tri2 = { rows[i][j], rows[i][j + 1], rows[i + 1][j + 1] };

                for (auto& tri : {tri1, tri2}) 
                {
                    P_triangolato.Cell2DsVertices.push_back(tri);
                    P_triangolato.Cell2DsId.push_back(nextFaceId++);

                    vector<unsigned int> edgeIds; //lista degli id dei lati del triangolo corrente
                    for (int k = 0; k < 3; k++) 

                    
                    { // itera su ogni vertice del triangolo corrente
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
            for (int k = 0; k < 3; k++) 
            {
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
    P_triangolato.Cell3DsVertices = {P_triangolato.Cell0DsId};
    P_triangolato.Cell3DsEdges = {P_triangolato.Cell1DsId};
    P_triangolato.Cell3DsFaces = {P_triangolato.Cell2DsId};

    return P_triangolato;
}


// ***************************************************************************
// Funzione per la traingolazione di classe 2
Polyhedron triangulateClass2(const Polyhedron& P, const unsigned int& b) 
{
    Polyhedron P1 = triangulateClass1(P,b);
    Polyhedron P2;

    // counter
    unsigned int nextPointId = 0;
    unsigned int nextEdgeId = 0;
    unsigned int nextFaceId = 0;

    // Mappa: chiave = id inizio e fine del vertice, valore = id lato
    map<std::pair<unsigned int, unsigned int>, unsigned int> edgeKeyToId; // mappa per controllare i duplicati dei lati
    map<tuple<double, double, double>, unsigned int> pointMap;   // mappa per controllare i duplicati dei punti

    // Copia i punti originali
    for (unsigned int pid = 0; pid < P1.NumCell0Ds; pid++) {
        Vector3d pt = P1.Cell0DsCoordinates.col(pid);
        P2.Cell0DsCoordinates.conservativeResize(3, nextPointId + 1);
        P2.Cell0DsCoordinates.col(nextPointId) = pt;
        P2.Cell0DsId.push_back(nextPointId);
        pointMap[{pt(0), pt(1), pt(2)}] = nextPointId;
        ++nextPointId;
    }

    // Baricentri delle facce
    vector<unsigned int> faceBarycenterIds(P1.NumCell2Ds);
    for (unsigned int fid = 0; fid < P1.NumCell2Ds; ++fid) 
    {
        const auto& verts = P1.Cell2DsVertices[fid];
        Vector3d barycenter = Vector3d::Zero();
        for (auto vid : verts) barycenter += P1.Cell0DsCoordinates.col(vid);
        barycenter /= verts.size();
    
        auto key = make_tuple(barycenter(0), barycenter(1), barycenter(2));

        if (!pointMap.count(key)) {
            P2.Cell0DsCoordinates.conservativeResize(3, nextPointId + 1);
            P2.Cell0DsCoordinates.col(nextPointId) = barycenter;
            P2.Cell0DsId.push_back(nextPointId);
            pointMap[key] = nextPointId;
            faceBarycenterIds[fid] = nextPointId++;
        } else {
            faceBarycenterIds[fid] = pointMap[key];
        }
    }


    // Mappa: chiave = id faccia di P, valore = vettore di id facce di P1 contenute
    std::map<unsigned int, std::vector<unsigned int>> facce_P_to_facce_P1;

    // Numero di triangoli generati per ogni faccia P nella triangolazione 1
    unsigned int num_tri_per_face = P1.Cell2DsVertices.size() / P.Cell2DsVertices.size();

    // Itera su tutte le facce del poliedro originale P
    for (unsigned int fid_P = 0; fid_P < P.Cell2DsVertices.size(); ++fid_P) 
    {
        unsigned int start_tri = fid_P * num_tri_per_face;
        unsigned int end_tri = start_tri + num_tri_per_face;

        for (unsigned int fid_P1 = start_tri; fid_P1 < end_tri; ++fid_P1) 
        {
            facce_P_to_facce_P1[fid_P].push_back(fid_P1);
        }
    }
   
    
    for (unsigned int eid = 0; eid < P1.Cell1DsId.size(); ++eid) {
        unsigned int v1 = P1.Cell1DsExtrema(0, eid);
        unsigned int v2 = P1.Cell1DsExtrema(1, eid);
        edgeKeyToId[std::minmax(v1, v2)] = P1.Cell1DsId[eid];
    }

    // Itera su tutti i triangoli di P1
    for (const auto& [fid_P, facceP1] : facce_P_to_facce_P1) 
    {
        // Mappa: chiave = id lato di P1, valore = vettore di id delle facce di P1 che lo contengono
        std::map<unsigned int, std::vector<unsigned int>> facce_per_lati;

        for (unsigned int fid_P1 : facceP1) {
            const std::vector<unsigned int>& edges = P1.Cell2DsEdges[fid_P1];

            for (unsigned int edge_id : edges) {
                auto it = facce_per_lati.find(edge_id);
                if (it != facce_per_lati.end()) {
                    std::vector<unsigned int>& faces_for_edge = it->second;
                    // aggiungo fid_P1 se non è già presente
                    if (std::find(faces_for_edge.begin(), faces_for_edge.end(), fid_P1) == faces_for_edge.end()) {
                        faces_for_edge.push_back(fid_P1);
                    }
                } else {
                    facce_per_lati[edge_id] = { fid_P1 };
                }
            }
        }

        // Collega ciascun baricentro ai vertici della faccia triangolare di classe 1
        for (unsigned int fid_P1 : facceP1) {
            unsigned int barycenter_id = faceBarycenterIds[fid_P1];
            const std::vector<unsigned int>& verts = P1.Cell2DsVertices[fid_P1];

            for (unsigned int vid : verts) {
                std::pair<unsigned int, unsigned int> edge_key = std::minmax(barycenter_id, vid);

                // Verifica se l'arco esiste già
                if (edgeKeyToId.find(edge_key) == edgeKeyToId.end()) {
                    unsigned int new_edge_id = P2.Cell1DsId.size();
                    P2.Cell1DsId.push_back(new_edge_id);
                    P2.Cell1DsExtrema.conservativeResize(2, new_edge_id + 1);
                    P2.Cell1DsExtrema(0, new_edge_id) = edge_key.first;
                    P2.Cell1DsExtrema(1, new_edge_id) = edge_key.second;
                    edgeKeyToId[edge_key] = new_edge_id;
                }
            }
        }        
        //distinzione tra caso in cui il lato appartiene ad una sola faccia e caso in cui il lato appartiene a due facce
        //caso in cui il lato appartiene ad una sola faccia
        //mappa temporanea: chiave = id faccia P1, valore = vector dei punti medi
        std::map<unsigned int, std::vector<unsigned int>> faceP1_to_midpoints;

        for (const auto& entry : facce_per_lati)
        {
            unsigned int edge_id = entry.first;
            const std::vector<unsigned int> faces = entry.second;

            // Caso: il lato appartiene a una sola faccia
            if (faces.size() == 1) {
                unsigned int fid_P1 = faces[0];
                unsigned int barycenter_id = faceBarycenterIds[fid_P1];

                unsigned int v1 = P1.Cell1DsExtrema(0, edge_id);
                unsigned int v2 = P1.Cell1DsExtrema(1, edge_id);

                // Calcolo punto medio
                Vector3d point1 = P1.Cell0DsCoordinates.col(v1);
                Vector3d point2 = P1.Cell0DsCoordinates.col(v2);
                Vector3d midpoint = 0.5 * (point1 + point2);

                // Controllo duplicati
                bool exists = false;
                unsigned int midpoint_id = 0;
                for (unsigned int pid = 0; pid < P2.Cell0DsCoordinates.cols(); ++pid) {
                    Vector3d existing = P2.Cell0DsCoordinates.col(pid);
                    if ((midpoint - existing).norm() < 1e-8) {
                        exists = true;
                        midpoint_id = pid;
                        break;
                    }
                }

                if (!exists) {
                    midpoint_id = P2.Cell0DsCoordinates.cols();
                    P2.Cell0DsCoordinates.conservativeResize(3, midpoint_id + 1);
                    P2.Cell0DsCoordinates.col(midpoint_id) = midpoint;
                    P2.Cell0DsId.push_back(midpoint_id);
                }

                faceP1_to_midpoints[fid_P1].push_back(midpoint_id);

                // Lati
                std::pair<unsigned int, unsigned int> edge_key1 = std::minmax(v1, midpoint_id);
                if (edgeKeyToId.find(edge_key1) == edgeKeyToId.end()) {
                    unsigned int new_edge_id = P2.Cell1DsId.size();
                    P2.Cell1DsId.push_back(new_edge_id);
                    P2.Cell1DsExtrema.conservativeResize(2, new_edge_id + 1);
                    P2.Cell1DsExtrema(0, new_edge_id) = edge_key1.first;
                    P2.Cell1DsExtrema(1, new_edge_id) = edge_key1.second;
                    edgeKeyToId[edge_key1] = new_edge_id;
                }

                std::pair<unsigned int, unsigned int> edge_key2 = std::minmax(v2, midpoint_id);
                if (edgeKeyToId.find(edge_key2) == edgeKeyToId.end()) {
                    unsigned int new_edge_id = P2.Cell1DsId.size();
                    P2.Cell1DsId.push_back(new_edge_id);
                    P2.Cell1DsExtrema.conservativeResize(2, new_edge_id + 1);
                    P2.Cell1DsExtrema(0, new_edge_id) = edge_key2.first;
                    P2.Cell1DsExtrema(1, new_edge_id) = edge_key2.second;
                    edgeKeyToId[edge_key2] = new_edge_id;
                }

                // Facce triangolari
                std::vector<unsigned int> triangle1 = { v1, midpoint_id, barycenter_id };
                std::vector<unsigned int> triangle2 = { v2, midpoint_id, barycenter_id };

                std::vector<unsigned int> edges1 = {
                    edgeKeyToId[std::minmax(v1, midpoint_id)],
                    edgeKeyToId[std::minmax(midpoint_id, barycenter_id)],
                    edgeKeyToId[std::minmax(barycenter_id, v1)]
                };
                std::vector<unsigned int> edges2 = {
                    edgeKeyToId[std::minmax(v2, midpoint_id)],
                    edgeKeyToId[std::minmax(midpoint_id, barycenter_id)],
                    edgeKeyToId[std::minmax(barycenter_id, v2)]
                };

                P2.Cell2DsId.push_back(P2.Cell2DsId.size());
                P2.Cell2DsVertices.push_back(triangle1);
                P2.Cell2DsEdges.push_back(edges1);

                P2.Cell2DsId.push_back(P2.Cell2DsId.size());
                P2.Cell2DsVertices.push_back(triangle2);
                P2.Cell2DsEdges.push_back(edges2);
            }


            // Caso: lato appartiene a due facce (da gestire se necessario)
            else if (faces.size() == 2) 
            {
                unsigned int fid1 = faces[0];
                unsigned int fid2 = faces[1];

                unsigned int bary1 = faceBarycenterIds[fid1];
                unsigned int bary2 = faceBarycenterIds[fid2];

                unsigned int v1 = P1.Cell1DsExtrema(0, edge_id);
                unsigned int v2 = P1.Cell1DsExtrema(1, edge_id);

                std::pair<unsigned int, unsigned int> edge_key = std::minmax(bary1, bary2);
                if (edgeKeyToId.find(edge_key) == edgeKeyToId.end()) {
                    unsigned int new_edge_id = P2.Cell1DsId.size();
                    P2.Cell1DsId.push_back(new_edge_id);
                    P2.Cell1DsExtrema.conservativeResize(2, new_edge_id + 1);
                    P2.Cell1DsExtrema(0, new_edge_id) = edge_key.first;
                    P2.Cell1DsExtrema(1, new_edge_id) = edge_key.second;
                    edgeKeyToId[edge_key] = new_edge_id;
                }

                // Due facce triangolari: (v1, bary1, bary2) e (v2, bary1, bary2)
                std::vector<unsigned int> triangle1 = { v1, bary1, bary2 };
                std::vector<unsigned int> triangle2 = { v2, bary1, bary2 };

                std::vector<unsigned int> edges1 = {
                    edgeKeyToId[std::minmax(v1, bary1)],
                    edgeKeyToId[std::minmax(v1, bary2)],
                    edgeKeyToId[std::minmax(bary1, bary2)]
                };
                std::vector<unsigned int> edges2 = {
                    edgeKeyToId[std::minmax(v2, bary1)],
                    edgeKeyToId[std::minmax(v2, bary2)],
                    edgeKeyToId[std::minmax(bary1, bary2)]
                };

                P2.Cell2DsId.push_back(P2.Cell2DsId.size());
                P2.Cell2DsVertices.push_back(triangle1);
                P2.Cell2DsEdges.push_back(edges1);

                P2.Cell2DsId.push_back(P2.Cell2DsId.size());
                P2.Cell2DsVertices.push_back(triangle2);
                P2.Cell2DsEdges.push_back(edges2);
            }
        }
    

        // Ora collega ciascun baricentro ai punti medi della stessa faccia
        for (const auto& entry : faceP1_to_midpoints) {
            unsigned int fid_P1 = entry.first;
            unsigned int barycenter_id = faceBarycenterIds[fid_P1];
            const std::vector<unsigned int>& midpoints = entry.second;

            for (unsigned int midpoint_id : midpoints) {
                // Controllo duplicati: cerca se il lato esiste già in P2
                bool edge_exists = false;
                for (unsigned int eid = 0; eid < P2.Cell1DsId.size(); ++eid) {
                    unsigned int e_v1 = P2.Cell1DsExtrema(0, eid);
                    unsigned int e_v2 = P2.Cell1DsExtrema(1, eid);
                    if ((e_v1 == barycenter_id && e_v2 == midpoint_id) ||
                        (e_v1 == midpoint_id && e_v2 == barycenter_id)) {
                        edge_exists = true;
                        break;
                    }
                }

                if (!edge_exists) {
                    unsigned int new_edge_id = P2.Cell1DsId.size();
                    P2.Cell1DsId.push_back(new_edge_id);
                    P2.Cell1DsExtrema.conservativeResize(2, new_edge_id + 1);
                    P2.Cell1DsExtrema(0, new_edge_id) = barycenter_id;
                    P2.Cell1DsExtrema(1, new_edge_id) = midpoint_id;
                }
            }
        }
    }

    //Aggiornamento strutture
    P2.NumCell0Ds = P2.Cell0DsId.size();
    P2.NumCell1Ds = P2.Cell1DsId.size();
    P2.NumCell2Ds = P2.Cell2DsId.size();

    P2.Cell3DsVertices.push_back(P2.Cell0DsId); 
    P2.Cell3DsEdges.push_back(P2.Cell1DsId); 
    P2.Cell3DsFaces.push_back(P2.Cell2DsId); 

    return P2;
}

// ***************************************************************************
// Funzione per la dualizzazione
Polyhedron Dualize(const Polyhedron& P_normale) 
{
    Polyhedron P_duale;
    // 1. Ogni faccia originale diventa un vertice del duale
    unsigned int num_faces = P_normale.Cell2DsVertices.size();
    P_duale.NumCell0Ds = num_faces;
    P_duale.Cell0DsId.resize(num_faces);
    P_duale.Cell0DsCoordinates.resize(3, num_faces);

    for (unsigned int i = 0; i < num_faces; i++)
    {
        P_duale.Cell0DsId[i] = i;

        const vector<unsigned int> face = P_normale.Cell2DsVertices[i];
        Vector3d centroid(0, 0, 0);
        for (unsigned int vid : face) 
        {
            centroid += P_normale.Cell0DsCoordinates.col(vid);
        }
        centroid = centroid / face.size();
        P_duale.Cell0DsCoordinates.col(i) = centroid;
    }

    // 2. Costruisci mappa edge -> facce (per trovare adiacenze tra facce)
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> edge_to_faces;
    for (unsigned int i = 0; i < num_faces; i++) 
    {
        const vector<unsigned int> face = P_normale.Cell2DsVertices[i];
        unsigned int n = face.size();
        for (unsigned int j = 0; j < n; j++) 
        {
            unsigned int v1 = face[j];
            unsigned int v2 = face[(j + 1) % n];
            auto edge = minmax(v1, v2);
            edge_to_faces[edge].push_back(i);
        }
    }

    // 3. Ogni coppia di facce adiacenti genera un lato nel duale
    unsigned int edge_id = 0;
    vector<Vector2i> temp_edges;

    for (const auto& [edge, faces] : edge_to_faces) 
    {
        if (faces.size() == 2) 
        {
            unsigned int f1 = faces[0];
            unsigned int f2 = faces[1];
            P_duale.Cell1DsId.push_back(edge_id++);
            temp_edges.push_back(Vector2i(f1, f2));
        }
    }

    P_duale.NumCell1Ds = edge_id;
    P_duale.Cell1DsExtrema.resize(2, edge_id);
    for (unsigned int i = 0; i < edge_id; ++i)
    {
        P_duale.Cell1DsExtrema.col(i) = temp_edges[i];
    }


    // 4. Ogni vertice originale diventa una faccia nel duale
    map<unsigned int, vector<unsigned int>> vertex_to_faces;
    for (unsigned int i = 0; i < num_faces; i++)
    {
        for (unsigned int vid : P_normale.Cell2DsVertices[i])
        {
            vertex_to_faces[vid].push_back(i);
        }
    }

    unsigned int face_id = 0;
    for (const auto& [vid, faces] : vertex_to_faces)
    {
        vector<unsigned int> face = faces;

        P_duale.Cell2DsId.push_back(face_id);
        P_duale.Cell2DsVertices.push_back(face);
        ++face_id;
    }
    P_duale.NumCell2Ds = face_id;

    // 5. Assegna l’unica cella 3D del duale
    P_duale.NumCell3Ds = 1;
    P_duale.Cell3DsFaces.resize(1);
    for (unsigned int i = 0; i < P_duale.NumCell2Ds; i++) 
    {
        P_duale.Cell3DsFaces[0].push_back(i);
    }

    // Aggiornamento Cell3D
    P_duale.Cell3DsVertices = {P_duale.Cell0DsId};
    P_duale.Cell3DsEdges = {P_duale.Cell1DsId};
    P_duale.Cell3DsFaces = {P_duale.Cell2DsId};

    // Riempimento Cell2DsEdges
    map<pair<unsigned int, unsigned int>, unsigned int> edge_lookup;
    for (unsigned int eid = 0; eid < P_duale.NumCell1Ds; eid++)
    {
        unsigned int v1 = P_duale.Cell1DsExtrema(0, eid);
        unsigned int v2 = P_duale.Cell1DsExtrema(1, eid);
        auto key = minmax(v1, v2);
        edge_lookup[key] = eid;
    }

    P_duale.Cell2DsEdges.resize(P_duale.NumCell2Ds);
    for (unsigned int fid = 0; fid < P_duale.NumCell2Ds; fid++) 
    {
        const vector<unsigned int> verts = P_duale.Cell2DsVertices[fid];
        vector<unsigned int> edge_ids;
        unsigned int n = verts.size();

        for (unsigned int i = 0; i < n; i++) 
        {
            unsigned int v1 = verts[i];
            unsigned int v2 = verts[(i + 1) % n];
            auto key = minmax(v1, v2);

            auto it = edge_lookup.find(key);
            if (it != edge_lookup.end()) 
            {
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
    return P_duale;   
}


// ***************************************************************************
Polyhedron projectPolyhedronOnSphere(const Polyhedron& P)
{
    Polyhedron projected = P;
    for (int i = 0; i < projected.Cell0DsCoordinates.cols(); i++)
    {
        projected.Cell0DsCoordinates.col(i).normalize(); 
    }
    return projected;
}


// ***************************************************************************
void ExportPolyhedron(const Polyhedron& P, const vector<Gedim::UCDProperty<double>>& points_properties, const vector<Gedim::UCDProperty<double>>& segments_properties)
{
    Gedim::UCDUtilities utilities;
    {
        utilities.ExportPoints("./Cell0Ds.inp",
                            P.Cell0DsCoordinates,
                            points_properties);
    }

    {
        utilities.ExportSegments("./Cell1Ds.inp",
                                P.Cell0DsCoordinates,
                                P.Cell1DsExtrema,
                                points_properties,
                                segments_properties);
    }
    Esporta_file(P);
}

// ***************************************************************************
void Esporta_file(const Polyhedron& P)
{
    //File: "Cell0Ds.txt"
    ofstream ofile1("Cell0Ds.txt");
    ofile1 << "Id;X;Y;Z\n"; 
    
    const MatrixXd scrivi_vertici = P.Cell0DsCoordinates;
    for(unsigned int id = 0; id < P.NumCell0Ds; id++)
        ofile1 << defaultfloat << id << ';' << scientific << setprecision(16) << scrivi_vertici(0,id) << ';' << scrivi_vertici(1,id) << ';' << scrivi_vertici(2,id) << '\n';
    ofile1.close();

    //Cell1Ds.txt
    ofstream ofile2("Cell1Ds.txt");
    ofile2 << "Id;Origin;End\n"; //header

    const MatrixXi scrivi_lati = P.Cell1DsExtrema;
    for(unsigned int id = 0; id < P.NumCell1Ds; id++)
        ofile2 << id << ';' << scrivi_lati(0,id) << ';' << scrivi_lati(1,id) << '\n';
    ofile2.close();

    //File: "Cell2Ds.txt"
    ofstream ofile3("Cell2Ds.txt");
    ofile3 << "Id;NumVertices;Vertices;NumEdges;Edges\n"; 
    
    const vector<vector<unsigned int>> scrivi_vertici_facce = P.Cell2DsVertices;
    const vector<vector<unsigned int>> scrivi_lati_facce = P.Cell2DsEdges;
    for(unsigned int id = 0; id < P.NumCell2Ds; id++)
    {
        ofile3 << id << ';' << scrivi_vertici_facce[id].size();
        for(unsigned int i = 0; i < scrivi_vertici_facce[id].size(); i++)
            ofile3 << ';' << scrivi_vertici_facce[id][i];
        
        ofile3 << ';' << scrivi_lati_facce[id].size();
        for(unsigned int j = 0; j < scrivi_lati_facce[id].size(); j++)
            ofile3 << ';' << scrivi_lati_facce[id][j];

        ofile3 << '\n';
    }
    ofile3.close();

    //File: "Cell3Ds.txt"
    ofstream ofile4("Cell3Ds.txt");
    ofile4 << "Id;NumVertices;Vertices;NumEdges;Edges;NumFaces;Faces\n";

    const vector<vector<unsigned int>> scrivi_vertici_poliedro = P.Cell3DsVertices;
    const vector<vector<unsigned int>> scrivi_lati_poliedro = P.Cell3DsEdges;
    const vector<vector<unsigned int>> scrivi_facce_poliedro = P.Cell3DsFaces;

    // visto che abbiamo un solo poliedro, riscriviamo direttamente i dati della struct senza il ciclo
    unsigned int id =  P.Cell3DsId;
    ofile4 << id << ';' << scrivi_vertici_poliedro[id].size();
    for(unsigned int i = 0; i < scrivi_vertici_poliedro[id].size(); i++)
        ofile4 << ';' << scrivi_vertici_poliedro[id][i];
    
    ofile4 << ';' << scrivi_lati_poliedro[id].size();
    for(unsigned int j = 0; j < scrivi_lati_poliedro[id].size(); j++)
        ofile4 << ';' << scrivi_lati_poliedro[id][j];

    ofile4 << ';' << scrivi_facce_poliedro[id].size();
    for(unsigned int j = 0; j < scrivi_facce_poliedro[id].size(); j++)
        ofile4 << ';' << scrivi_facce_poliedro[id][j];

    ofile4 << '\n';

    ofile4.close();

}

// ***************************************************************************
vector<unsigned int> Cammini_minimi(const Polyhedron& P, const unsigned int& v1, const unsigned int& v2)
{
    const MatrixXi edges = P.Cell1DsExtrema;
    const MatrixXd coords = P.Cell0DsCoordinates;

    // Costruisco lista adiacenza con pesi 
    vector<vector<pair<unsigned int, double>>> lista_adiacenza(P.NumCell0Ds); 
    for(unsigned int j = 0; j < edges.cols(); j++)
    {
        unsigned int inizio = edges(0, j);
        unsigned int fine = edges(1, j);
        Vector3d p1 = coords.col(inizio);
        Vector3d p2 = coords.col(fine);
        double distanza = (p1 - p2).norm();

        lista_adiacenza[inizio].emplace_back(fine, distanza);
        lista_adiacenza[fine].emplace_back(inizio, distanza);
    }

    // Implementazione di Dijkstra
    const unsigned int N = P.NumCell0Ds;
    vector<double> distanza(N, numeric_limits<double>::max());
    vector<unsigned int> parent(N, N); 
    priority_queue<pair<double, unsigned int>, vector<pair<double, unsigned int>>, greater<pair<double, unsigned int>>> heap;

    distanza[v1] = 0.0;
    heap.emplace(0.0, v1);

    while(!heap.empty())
    {
        auto [dist_corrente, nodo_corrente] = heap.top();
        heap.pop();

        if(nodo_corrente == v2) 
            break; // trovato cammino più breve fino a v2

        for(auto [vicino, peso] : lista_adiacenza[nodo_corrente])
        {
            double nuova_distanza = dist_corrente + peso;
            if(nuova_distanza < distanza[vicino])
            {
                distanza[vicino] = nuova_distanza;
                parent[vicino] = nodo_corrente;
                heap.emplace(nuova_distanza, vicino);
            }
        }
    }

    // Ricostruzione percorso
    vector<unsigned int> path;
    if(distanza[v2] == numeric_limits<double>::max())
    {
        cout << "Nessun percorso trovato tra " << v1 << " e " << v2 << endl;
        return path;
    }

    for(unsigned int node = v2; node != N; node = parent[node])
    {
        path.push_back(node);
    }
    reverse(path.begin(), path.end());

    cout<<"Il cammino minimo che collega "<< v1 <<" e "<< v2 <<" ha "<<(path.size()-1)<< " lati "<< endl;
    cout << "Il cammino minimo che collega " << v1 << " e " << v2 << " ha lunghezza: " << distanza[v2] << endl;

    return path;
}

}

