#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace PolyhedronLibrary
{
bool ImportMesh(Polyhedron& mesh)
{

    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

    return true;

}
// ***************************************************************************
// importo le celle 0D
bool ImportCell0Ds(Polyhedron& mesh)
{
    ifstream file("./Cell0Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);
        string token;

        unsigned int id;
        Vector2d coord;

        // Adesso leggi i dati separati da ';'
        getline(converter, token, ';');
        id = static_cast<unsigned int>(stoi(token));

        getline(converter, token, ';');

        getline(converter, token, ';');
        coord(0) = stod(token);

        getline(converter, token, ';');
        coord(1) = stod(token);

        mesh.Cell0DsId.push_back(id);

        mesh.Cell0DsCoordinates(0, id) = coord(0);
        mesh.Cell0DsCoordinates(1, id) = coord(1);

        //new
        const auto insert = mesh.IdCell0Ds.find(id);
            if(insert == mesh.IdCell0Ds.end())
            {
                mesh.IdCell0Ds.insert({id, {coord(0), coord(1)}});
            }
            else
            {
                insert->second.push_back(id);
            }

    }

    return true;
}
// ***************************************************************************
//Qui viene fatta la stessa cosa con le celle 1D, con 1 id, 1 marker, 1 origine e 1 fine (del segmento)
bool ImportCell1Ds(Polyhedron& mesh)
{
    ifstream file("./Cell1Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds); // 2 righe, NumCell1Ds colonne

    for (const string& line : listLines)
    {
        istringstream converter(line);
        string token;

        unsigned int id;
        unsigned int vertex0;
        unsigned int vertex1;

        getline(converter, token, ';');
        id = static_cast<unsigned int>(stoi(token));

        getline(converter, token, ';');

        getline(converter, token, ';');
        vertex0 = static_cast<unsigned int>(stoi(token));

        getline(converter, token, ';');
        vertex1 = static_cast<unsigned int>(stoi(token));

        mesh.Cell1DsExtrema(0, id) = vertex0;
        mesh.Cell1DsExtrema(1, id) = vertex1;

        mesh.Cell1DsId.push_back(id);

    }

    return true;
}
// ***************************************************************************
//Qui viene fatta la stessa cosa con le celle 2D, con 1 id, 3 vertici e 3 lati
bool ImportCell2Ds(Polyhedron& mesh)
{
    ifstream file;
    file.open("./Cell2Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines; // creo una lista che conterrà tutte le righe del file salvate come stringhe
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close(); // non mi serve più il file quindi lo chiudo

    // tolgo l'intestazione
    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size(); //conta tutte le righe rimanenti, che corrsipondono al numero di POLIGONI (celle 2D)

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds); //spazio per id dei triangoli
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds); //spazio per vertici dei triangoli
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds); //spazio per lati dei triangoli

    for (const string& line : listLines)
    {
        istringstream converter(line); 

        // Leggo i primi tre campi della riga
        unsigned int id, marker, numVertices;
        char delimiter;

        converter >> id >> delimiter >> marker >> delimiter >> numVertices;

        vector<unsigned int> vertices(numVertices); //creo un vettore che contiene i vertici (che sono di un numero pari a numVertices)
        for(unsigned int i = 0; i < numVertices; i++)
            converter >> delimiter >> vertices[i];

        unsigned int numEdges; //leggo il numero di lati
        converter >> delimiter >> numEdges;
        //std::cout << line << std::endl;
        vector<unsigned int> edges(numEdges); //creo un vettore che contiene i lati (che sono di un numero pari a numEdges)
        for(unsigned int i = 0; i < numEdges; i++)
            converter >> delimiter >> edges[i];

        
        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);
    }

    return true;
}

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
}

// ***************************************************************************
// Funzione per la triangolazione
