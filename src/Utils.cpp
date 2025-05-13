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

}