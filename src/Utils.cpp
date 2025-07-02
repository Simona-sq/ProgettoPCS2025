#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <set>
#include "Utils.hpp"
#include "UCDUtilities.hpp"
#include "Polyhedron.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedronLibrary 
{
    // ***************************************************************************
    // FUNZIONE SalvataggioDati
    // Costruisco i dati geometrici (vertici e facce) di un solido platonico
    pair<vector<Vector3d>, vector<vector<unsigned int>>> SalvataggioDati(unsigned int& q) 
    {   
        // Inizializzo due vettori vuoti per salvare i dati
        vector<Vector3d> vertici;
        vector<vector<unsigned int>> facce;

        // Controllo il valore di q per determinare quale solido costruire
        if (q == 3) // Tetraedro 
        {  
            vertici = {
                {0.0, 0.0, sqrt(3)},
                {sqrt(8.0/3.0), 0.0, -1.0/sqrt(3)},
                {-sqrt(2.0/3.0), sqrt(2.0), -1.0/sqrt(3)},
                {-sqrt(2.0/3.0), -sqrt(2.0), -1.0/sqrt(3)}
            };
            for (auto& v : vertici) v = v.normalized();
            facce = {
                {0, 1, 2},
                {0, 2, 3},
                {0, 3, 1},
                {1, 3, 2}
            };

        } 
        else if (q == 4) // Ottaedro
        { 
            vertici = {
                {1, 0, 0}, {-1, 0, 0},
                {0, 1, 0}, {0, -1, 0},
                {0, 0, 1}, {0, 0, -1}
            };
            for (auto& v : vertici) v = v.normalized();
            facce = {
                {0, 4, 2}, {2, 4, 1}, {1, 4, 3}, {3, 4, 0},
                {0, 2, 5}, {2, 1, 5}, {1, 3, 5}, {3, 0, 5}
            };
        } 
        else if (q == 5)  // Icosaedro
        {
            vector<Vector3d> coordinate = {
                { 0.0, 0.0, 1.175571}, { 1.051462, 0.0, 0.5257311}, {0.3249197, 1.0, 0.5257311},
                { -0.8506508, 0.618034, 0.5257311}, {-0.8506508, -0.618034, 0.5257311}, { 0.3249197, -1.0, 0.5257311}, 
                { 0.8506508, 0.618034, -0.5257311}, { 0.8506508, -0.618034, -0.5257311}, {-0.3249197, 1.0, -0.5257311}, 
                {-1.051462, 0.0, -0.5257311}, {-0.3249197, -1.0, -0.5257311}, {0.0, 0.0, -1.175571}
            };
            for (auto& v : coordinate) vertici.push_back(v.normalized());
            facce = {
                {0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 5, 1},
                {1, 5, 7}, {1, 7, 6}, {1, 6, 2}, {2, 6, 8}, {2, 8, 3},
                {3, 8, 9}, {3, 9, 4}, {4, 9, 10}, {4, 10, 5}, {5, 10, 7},
                {6, 7, 11}, {6, 11, 8}, {7, 10, 11}, {8, 11, 9}, {9, 11, 10}
            };
        }

        return {vertici, facce};
    }


    // ***************************************************************************
    // FUNZIONE CreazioneSolidoPlatonico
    // Costruisco e restituisco l'oggetto Polyhedron che rappresenta un solido platonico
    Polyhedron CreazioneSolidoPlatonico(unsigned int& q) 
    {
        Polyhedron P;
        // Recupero dati geometrici
        auto [vertici, facce] = SalvataggioDati(q);

        // Riempio le info sui vertici (Cell0Ds)
        P.NumCell0Ds = vertici.size();
        P.Cell0DsCoordinates = MatrixXd(3, P.NumCell0Ds);
        for (unsigned int i = 0; i < vertici.size(); i++) 
        {
            P.Cell0DsId.push_back(i);
            P.Cell0DsCoordinates.col(i) = vertici[i];
        }
        
        // Riempio le info sulle facce (Cell2Ds)
        // Creo lati univoci con orientamento coerente
        map<pair<unsigned int, unsigned int>, unsigned int> MappaLati;
        unsigned int contatore_lati = 0;

        // Per ogni faccia: registro l’ID della faccia, salvo i vertici della faccia ed estraggo tutti i lati, tenendo conto dell’orientamento
        for (unsigned int fid = 0; fid < facce.size(); fid++) 
        {
            const auto& f = facce[fid];
            P.Cell2DsId.push_back(fid);
            P.Cell2DsVertices.push_back(f);

            vector<unsigned int> IdLati;
            for (unsigned int i = 0; i < f.size(); i++)
            {
                unsigned int v1 = f[i];
                unsigned int v2 = f[(i + 1) % f.size()];
                pair<unsigned int, unsigned int> key = {v1, v2};
                pair<unsigned int, unsigned int> rev_key = {v2, v1};

                if (MappaLati.find(key) == MappaLati.end() && MappaLati.find(rev_key) == MappaLati.end()) 
                {
                    MappaLati[key] = contatore_lati;
                    P.Cell1DsId.push_back(contatore_lati);
                    contatore_lati++;
                }

                // Se il lato non è stato registrato (né diretto né inverso), lo aggiungo alla mappa con un ID
                if (MappaLati.find(key) != MappaLati.end()) 
                {
                    IdLati.push_back(MappaLati[key]);
                } 
                else 
                {
                    IdLati.push_back(MappaLati[rev_key]);
                }
            }
            P.Cell2DsEdges.push_back(IdLati);
        }

        // Riempio le info sui lati (Cell1Ds)
        P.NumCell1Ds = MappaLati.size();
        // Alloco una matrice 2×N con estremi di ciascun lato e la riempio
        P.Cell1DsExtrema = MatrixXi(2, P.NumCell1Ds);
        for (const auto& [key, eid] : MappaLati) 
        {
            P.Cell1DsExtrema(0, eid) = key.first;
            P.Cell1DsExtrema(1, eid) = key.second;
        }

        // Riempio le info sul poliedro (Cell3Ds)
        P.NumCell2Ds = facce.size();
        P.Cell3DsVertices = P.Cell0DsId;
        P.Cell3DsEdges = P.Cell1DsId;
        P.Cell3DsFaces = P.Cell2DsId;

        return P;
    }


    // ***************************************************************************
    // FUNZIONE TriangolazioneClasse1
    // Creo la triangolazione di classe 1 del solido platonico di partenza
    Polyhedron TriangolazioneClasse1(const Polyhedron& P, const unsigned int& t_value)
    {
        // Inizializzazione varibaili e strutture
        Polyhedron P1;
        
        unsigned int id_punto_successivo = 0;
        unsigned int id_lato_successivo = 0;
        unsigned int id_faccia_successiva = 0;

        map<pair<unsigned int, unsigned int>, unsigned int> MappaLati; // mappa che tiene traccia dei lati già esistenti evitando duplicati
        map<tuple<double, double, double>, unsigned int> MappaPunti; // mappa che tiene traccia dei punti già esistenti evitando duplicati

        // Itero su ogni faccia del solido platonico di partenza
        for (size_t fid = 0; fid < P.Cell2DsVertices.size(); fid++) 
        {
            const auto& faccia = P.Cell2DsVertices[fid];

            // Salvo gli id dei vertici A B C del solido platonico di partenza
            unsigned int A = faccia[0]; 
            unsigned int B = faccia[1];
            unsigned int C = faccia[2];

            // Salvo le coordinate dei vertici A B C del solido platonico di partenza
            Vector3d vA = P.Cell0DsCoordinates.col(A);
            Vector3d vB = P.Cell0DsCoordinates.col(B);
            Vector3d vC = P.Cell0DsCoordinates.col(C);

            vector<vector<unsigned int>> rows;

            // Costruisco i punti interni con interpolazione baricentrica
            for (unsigned int i = 0; i < t_value+1; i++) // ogni iterazione rappresenta una riga di punti nel triangolo
            { 
                vector<unsigned int> row; // contiene gli id dei punti generati in quella riga

                for (unsigned int j = 0; j <= i; j++) // ogni iterazione rappresenta un punto sulla riga
                { 
                    // Trovo i parametri per dividere il lato i t_value parti
                    double alpha = double(t_value - i) / t_value; 
                    double beta  = double(i - j) / t_value;
                    double gamma = double(j) / t_value;
                    Vector3d point = alpha * vA + beta * vB + gamma * vC; // coordinate del nuovo punto
                    
                    // Verifico i duplicati 
                    auto key = make_tuple(point(0), point(1), point(2));
                    auto it = MappaPunti.find(key);
                    if (it != MappaPunti.end()) // caso 1: il punto esiste
                    {
                        row.push_back(it->second);
                    } 
                    else // caso 2: il punto è nuovo
                    {
                        // Salvo il nuovo punto
                        P1.Cell0DsCoordinates.conservativeResize(3, id_punto_successivo + 1);
                        P1.Cell0DsCoordinates.col(id_punto_successivo) = point;
                        P1.Cell0DsId.push_back(id_punto_successivo);
                        row.push_back(id_punto_successivo);
                        MappaPunti[key] = id_punto_successivo;
                        ++id_punto_successivo;
                    }  
                }
                rows.push_back(row);
            }

            // Creo i triangoli
            for (unsigned int i = 0; i < t_value; i++)  // ogni iterazione rappresenta una riga di triangoli che verranno creati
            { 
                // Parte 1: Triangolazione di tutta la faccia tranne il triangolo in cima, a partire dalla linea con due punti
                for (unsigned int j = 0; j < i; j++)  // ogni iterazione rappresenta un triangolo che verrà creato nella riga corrente
                { 
                    vector<unsigned int> triangolo1 = { rows[i][j], rows[i + 1][j], rows[i + 1][j + 1] };
                    vector<unsigned int> triangolo2 = { rows[i][j], rows[i][j + 1], rows[i + 1][j + 1] };

                    for (auto& triangolo : {triangolo1, triangolo2}) 
                    {
                        P1.Cell2DsVertices.push_back(triangolo);
                        P1.Cell2DsId.push_back(id_faccia_successiva++);

                        vector<unsigned int> IdLati; // lista degli id dei lati del triangolo corrente

                        for (int k = 0; k < 3; k++) // itera su ogni vertice del triangolo corrente
                        { 
                            auto v1 = triangolo[k]; // vertice corrente
                            auto v2 = triangolo[(k + 1) % 3]; // vertice successivo + condizione per la chiusura del tiangolo
                            auto key = minmax(v1, v2); // ordina la coppia per evitare duplicati
                            if (MappaLati.find(key) == MappaLati.end()) 
                            {
                                MappaLati[key] = id_lato_successivo++;
                                P1.Cell1DsId.push_back(MappaLati[key]); // aggiunge l'id del nuovo lato a Cell1DsId
                            }
                            IdLati.push_back(MappaLati[key]);
                        }
                        P1.Cell2DsEdges.push_back(IdLati);
                    }
                }

                // Parte 2: Triangolazione del triangolo in cima
                vector<unsigned int> triangolo = { rows[i][i], rows[i + 1][i], rows[i + 1][i + 1] };
                P1.Cell2DsVertices.push_back(triangolo);
                P1.Cell2DsId.push_back(id_faccia_successiva++);

                vector<unsigned int> IdLati;
                for (int k = 0; k < 3; k++) 
                {
                    auto v1 = triangolo[k];
                    auto v2 = triangolo[(k + 1) % 3];
                    auto key = minmax(v1, v2);
                    if (MappaLati.find(key) == MappaLati.end()) 
                    {
                        MappaLati[key] = id_lato_successivo++;
                        P1.Cell1DsId.push_back(MappaLati[key]);
                    }
                    IdLati.push_back(MappaLati[key]);
                }
                P1.Cell2DsEdges.push_back(IdLati);
            }
        }

        // Aggiorno gli attributi del poliedro
        P1.NumCell0Ds = id_punto_successivo;
        P1.NumCell1Ds = MappaLati.size();
        P1.NumCell2Ds = id_faccia_successiva;

        // la matrice Cell1DsExtrema viene riempita con i valori aggiornati di edgeMap
        P1.Cell1DsExtrema = MatrixXi(2, P1.NumCell1Ds); 
        for (const auto& [key, eid] : MappaLati) 
        {
            P1.Cell1DsExtrema(0, eid) = key.first;
            P1.Cell1DsExtrema(1, eid) = key.second;
        }

        P1.Cell3DsVertices = P1.Cell0DsId;
        P1.Cell3DsEdges = P1.Cell1DsId;
        P1.Cell3DsFaces = P1.Cell2DsId;

        return P1;
    }


    // ***************************************************************************
    // FUNZIONE TriangolazioneClasse2
    // Creo la triangolazione di classe 2 del solido platonico di partenza
    Polyhedron TriangolazioneClasse2(const Polyhedron& P, const unsigned int& b) 
    {   
        // Applico la triangolazione di classe 1 per iniziare la costruzione
        Polyhedron P1 = TriangolazioneClasse1(P,b);
        Polyhedron P2;

        // Contatore dei nuovi punti
        unsigned int id_nuovo_punto = 0;
 
        map<pair<unsigned int, unsigned int>, unsigned int> MappaLati; // mappa per controllare i duplicati dei lati: chiave = id inizio e fine del vertice, valore = id lato
        map<tuple<double, double, double>, unsigned int> MappaPunti;   // mappa per controllare i duplicati dei punti

        // Copio i punti originali da P1 in P2
        for (unsigned int pid = 0; pid < P1.NumCell0Ds; pid++) 
        {
            Vector3d punto = P1.Cell0DsCoordinates.col(pid);
            P2.Cell0DsCoordinates.conservativeResize(3, id_nuovo_punto + 1);
            P2.Cell0DsCoordinates.col(id_nuovo_punto) = punto;
            P2.Cell0DsId.push_back(id_nuovo_punto);
            MappaPunti[{punto(0), punto(1), punto(2)}] = id_nuovo_punto;
            ++id_nuovo_punto;
        }

        // Calcolo dei baricentri delle facce di P1
        vector<unsigned int> IdBaricentri(P1.NumCell2Ds);
        for (unsigned int fid = 0; fid < P1.NumCell2Ds; ++fid) 
        {
            const auto& vertici = P1.Cell2DsVertices[fid];
            Vector3d baricentro = Vector3d::Zero();
            for (auto vid : vertici) baricentro += P1.Cell0DsCoordinates.col(vid);
            baricentro /= vertici.size();
        
            auto key = make_tuple(baricentro(0), baricentro(1), baricentro(2));

            if (!MappaPunti.count(key)) 
            {
                P2.Cell0DsCoordinates.conservativeResize(3, id_nuovo_punto + 1);
                P2.Cell0DsCoordinates.col(id_nuovo_punto) = baricentro;
                P2.Cell0DsId.push_back(id_nuovo_punto);
                MappaPunti[key] = id_nuovo_punto;
                IdBaricentri[fid] = id_nuovo_punto++;
            }
            else
            {
                IdBaricentri[fid] = MappaPunti[key];
            }
        }

        // Mappa: chiave = id faccia di P, valore = vettore di id facce di P1 contenute
        map<unsigned int, vector<unsigned int>> FacceP1_inP;

        // Numero di triangoli generati per ogni faccia P nella triangolazione 1
        unsigned int num_tri_per_face = P1.Cell2DsVertices.size() / P.Cell2DsVertices.size();

        for (unsigned int fid_P = 0; fid_P < P.Cell2DsVertices.size(); fid_P++)  // itero su tutte le facce del poliedro originale P
        {
            unsigned int start_tri = fid_P * num_tri_per_face;
            unsigned int end_tri = start_tri + num_tri_per_face;

            for (unsigned int fid_P1 = start_tri; fid_P1 < end_tri; fid_P1++) 
            {
                FacceP1_inP[fid_P].push_back(fid_P1);
            }
        }
    
        // Costruzione mappa lati da P1
        for (unsigned int eid = 0; eid < P1.Cell1DsId.size(); eid++) 
        {
            unsigned int v1 = P1.Cell1DsExtrema(0, eid);
            unsigned int v2 = P1.Cell1DsExtrema(1, eid);
            MappaLati[minmax(v1, v2)] = P1.Cell1DsId[eid];
        }

        // Itera su tutti i triangoli di P1
        for (const auto& [fid_P, facceP1] : FacceP1_inP) 
        {
            // Mappa: chiave = id lato di P1, valore = vettore di id delle facce di P1 che lo contengono
            map<unsigned int, vector<unsigned int>> FaccePerLati;

            for (unsigned int fid_P1 : facceP1)
            {
                const vector<unsigned int>& lati = P1.Cell2DsEdges[fid_P1];

                for (unsigned int id_lato : lati) {
                    auto it = FaccePerLati.find(id_lato);
                    if (it != FaccePerLati.end()) 
                    {
                        vector<unsigned int>& facce_per_lati = it->second;
                        // aggiungo fid_P1 se non è già presente
                        if (find(facce_per_lati.begin(), facce_per_lati.end(), fid_P1) == facce_per_lati.end()) {
                            facce_per_lati.push_back(fid_P1);
                        }
                    } 
                    else 
                    {
                        FaccePerLati[id_lato] = { fid_P1 };
                    }
                }
            }

            // Collega ciascun baricentro ai vertici della faccia triangolare di classe 1
            for (unsigned int fid_P1 : facceP1) {
                unsigned int id_baricentro = IdBaricentri[fid_P1];
                const vector<unsigned int>& vertici = P1.Cell2DsVertices[fid_P1];

                for (unsigned int vid : vertici)
                {
                    pair<unsigned int, unsigned int> key_lato = minmax(id_baricentro, vid);

                    // Verifica se l'arco esiste già
                    if (MappaLati.find(key_lato) == MappaLati.end())
                    {
                        unsigned int id_nuovo_lato = P2.Cell1DsId.size();
                        P2.Cell1DsId.push_back(id_nuovo_lato);
                        P2.Cell1DsExtrema.conservativeResize(2, id_nuovo_lato + 1);
                        P2.Cell1DsExtrema(0, id_nuovo_lato) = key_lato.first;
                        P2.Cell1DsExtrema(1, id_nuovo_lato) = key_lato.second;
                        MappaLati[key_lato] = id_nuovo_lato;
                    }
                }
            }     

            // Gestione punti medi e costruzione facce triangolari
            map<unsigned int, vector<unsigned int>> PuntiMediP1;  // mappa temporanea: chiave = id faccia P1, valore = vector dei punti medi

            for (const auto& riga : FaccePerLati)
            {
                unsigned int id_lato = riga.first;
                const vector<unsigned int>& facce = riga.second;

                // Parte 1: il lato appartiene ad una sola faccia (sul bordo)
                if (facce.size() == 1) {
                    unsigned int fid_P1 = facce[0];
                    unsigned int id_baricentro = IdBaricentri[fid_P1];

                    unsigned int v1 = P1.Cell1DsExtrema(0, id_lato);
                    unsigned int v2 = P1.Cell1DsExtrema(1, id_lato);

                    // Calcolo punto medio
                    Vector3d punto1 = P1.Cell0DsCoordinates.col(v1);
                    Vector3d punto2 = P1.Cell0DsCoordinates.col(v2);
                    Vector3d punto_medio = 0.5 * (punto1 + punto2);

                    // Controllo duplicati
                    bool esiste = false;
                    unsigned int id_punto_medio = 0;
                    for (unsigned int pid = 0; pid < P2.Cell0DsCoordinates.cols(); pid++)
                    {
                        Vector3d esistente = P2.Cell0DsCoordinates.col(pid);
                        if ((punto_medio - esistente).norm() < 2.2e-16)
                        {
                            esiste = true;
                            id_punto_medio = pid;
                            break;
                        }
                    }

                    if (!esiste)
                    {
                        id_punto_medio = P2.Cell0DsCoordinates.cols();
                        P2.Cell0DsCoordinates.conservativeResize(3, id_punto_medio + 1);
                        P2.Cell0DsCoordinates.col(id_punto_medio) = punto_medio;
                        P2.Cell0DsId.push_back(id_punto_medio);
                    }

                    PuntiMediP1[fid_P1].push_back(id_punto_medio);

                    // Aggiungo i lati
                    pair<unsigned int, unsigned int> key1_lato = minmax(v1, id_punto_medio);
                    if (MappaLati.find(key1_lato) == MappaLati.end())
                    {
                        unsigned int id_nuovo_lato = P2.Cell1DsId.size();
                        P2.Cell1DsId.push_back(id_nuovo_lato);
                        P2.Cell1DsExtrema.conservativeResize(2, id_nuovo_lato + 1);
                        P2.Cell1DsExtrema(0, id_nuovo_lato) = key1_lato.first;
                        P2.Cell1DsExtrema(1, id_nuovo_lato) = key1_lato.second;
                        MappaLati[key1_lato] = id_nuovo_lato;
                    }

                    pair<unsigned int, unsigned int> key2_lato = minmax(v2, id_punto_medio);
                    if (MappaLati.find(key2_lato) == MappaLati.end())
                    {
                        unsigned int id_nuovo_lato = P2.Cell1DsId.size();
                        P2.Cell1DsId.push_back(id_nuovo_lato);
                        P2.Cell1DsExtrema.conservativeResize(2, id_nuovo_lato + 1);
                        P2.Cell1DsExtrema(0, id_nuovo_lato) = key2_lato.first;
                        P2.Cell1DsExtrema(1, id_nuovo_lato) = key2_lato.second;
                        MappaLati[key2_lato] = id_nuovo_lato;
                    }

                    // Creo due nuove facce triangolari
                    vector<unsigned int> triangolo1 = { v1, id_punto_medio, id_baricentro };
                    vector<unsigned int> triangolo2 = { v2, id_punto_medio, id_baricentro };

                    vector<unsigned int> lati1 = {
                        MappaLati[minmax(v1, id_punto_medio)],
                        MappaLati[minmax(id_punto_medio, id_baricentro)],
                        MappaLati[minmax(id_baricentro, v1)]
                    };
                    vector<unsigned int> lati2 = {
                        MappaLati[minmax(v2, id_punto_medio)],
                        MappaLati[minmax(id_punto_medio, id_baricentro)],
                        MappaLati[minmax(id_baricentro, v2)]
                    };

                    P2.Cell2DsId.push_back(P2.Cell2DsId.size());
                    P2.Cell2DsVertices.push_back(triangolo1);
                    P2.Cell2DsEdges.push_back(lati1);

                    P2.Cell2DsId.push_back(P2.Cell2DsId.size());
                    P2.Cell2DsVertices.push_back(triangolo2);
                    P2.Cell2DsEdges.push_back(lati2);
                }


                // Parte 2: il lato appartiene a due facce (interno)
                else if (facce.size() == 2) 
                {
                    unsigned int fid1 = facce[0];
                    unsigned int fid2 = facce[1];

                    unsigned int baricentro1 = IdBaricentri[fid1];
                    unsigned int baricentro2 = IdBaricentri[fid2];

                    unsigned int v1 = P1.Cell1DsExtrema(0, id_lato);
                    unsigned int v2 = P1.Cell1DsExtrema(1, id_lato);

                    pair<unsigned int, unsigned int> key_lato = minmax(baricentro1, baricentro2);
                    if (MappaLati.find(key_lato) == MappaLati.end())
                    {
                        unsigned int id_nuovo_lato = P2.Cell1DsId.size();
                        P2.Cell1DsId.push_back(id_nuovo_lato);
                        P2.Cell1DsExtrema.conservativeResize(2, id_nuovo_lato + 1);
                        P2.Cell1DsExtrema(0, id_nuovo_lato) = key_lato.first;
                        P2.Cell1DsExtrema(1, id_nuovo_lato) = key_lato.second;
                        MappaLati[key_lato] = id_nuovo_lato;
                    }

                    // Creo due facce triangolari: (v1, bary1, bary2) e (v2, bary1, bary2)
                    vector<unsigned int> triangolo1 = { v1, baricentro1, baricentro2 };
                    vector<unsigned int> triangolo2 = { v2, baricentro1, baricentro2 };

                    vector<unsigned int> lati1 = {
                        MappaLati[minmax(v1, baricentro1)],
                        MappaLati[minmax(v1, baricentro2)],
                        MappaLati[minmax(baricentro1, baricentro2)]
                    };
                    vector<unsigned int> lati2 = {
                        MappaLati[minmax(v2, baricentro1)],
                        MappaLati[minmax(v2, baricentro2)],
                        MappaLati[minmax(baricentro1, baricentro2)]
                    };

                    P2.Cell2DsId.push_back(P2.Cell2DsId.size());
                    P2.Cell2DsVertices.push_back(triangolo1);
                    P2.Cell2DsEdges.push_back(lati1);

                    P2.Cell2DsId.push_back(P2.Cell2DsId.size());
                    P2.Cell2DsVertices.push_back(triangolo2);
                    P2.Cell2DsEdges.push_back(lati2);
                }
            }

            // Collego ciascun baricentro ai punti medi della stessa faccia
            for (const auto& riga : PuntiMediP1)
            {
                unsigned int fid_P1 = riga.first;
                unsigned int id_baricentro = IdBaricentri[fid_P1];
                const vector<unsigned int>& punti_medi = riga.second;

                for (unsigned int id_punto_medio : punti_medi)
                {
                    // Controllo duplicati: cerco se il lato esiste già in P2
                    bool lato_esiste = false;
                    for (unsigned int eid = 0; eid < P2.Cell1DsId.size(); eid++)
                    {
                        unsigned int origine_lato = P2.Cell1DsExtrema(0, eid);
                        unsigned int fine_lato = P2.Cell1DsExtrema(1, eid);
                        if ((origine_lato == id_baricentro && fine_lato == id_punto_medio) ||
                            (origine_lato == id_punto_medio && fine_lato == id_baricentro))
                            {
                            lato_esiste = true;
                            break;
                        }
                    }

                    if (!lato_esiste)
                    {
                        unsigned int id_nuovo_lato = P2.Cell1DsId.size();
                        P2.Cell1DsId.push_back(id_nuovo_lato);
                        P2.Cell1DsExtrema.conservativeResize(2, id_nuovo_lato + 1);
                        P2.Cell1DsExtrema(0, id_nuovo_lato) = id_baricentro;
                        P2.Cell1DsExtrema(1, id_nuovo_lato) = id_punto_medio;
                    }
                }
            }
        }

        // Aggiorno struct
        P2.NumCell0Ds = P2.Cell0DsId.size();
        P2.NumCell1Ds = P2.Cell1DsId.size();
        P2.NumCell2Ds = P2.Cell2DsId.size();

        P2.Cell3DsVertices = P2.Cell0DsId; 
        P2.Cell3DsEdges = P2.Cell1DsId; 
        P2.Cell3DsFaces = P2.Cell2DsId; 

        return P2;
    }


    // ***************************************************************************
    // FUNZIONE Dualizzazione
    // Creo la dualizzazione
    Polyhedron Dualizzazione(const Polyhedron& P) 
    {
        Polyhedron P_duale;

        unsigned int numero_facce = P.Cell2DsVertices.size();
        P_duale.NumCell0Ds = numero_facce;
        P_duale.Cell0DsId.resize(numero_facce);
        P_duale.Cell0DsCoordinates.resize(3, numero_facce);

        // Costruisco i vertici del duale (baricentri delle facce di P)
        for (unsigned int i = 0; i < numero_facce; i++)
        {
            P_duale.Cell0DsId[i] = i;
            const vector<unsigned int>& faccia = P.Cell2DsVertices[i];
            Vector3d baricentro(0, 0, 0);
            for (unsigned int vid : faccia) 
            {
                baricentro += P.Cell0DsCoordinates.col(vid);
            }
            baricentro = baricentro / faccia.size();
            P_duale.Cell0DsCoordinates.col(i) = baricentro;
        }

        // Costruisco la mappa: lati -> facce (per trovare adiacenze tra facce)
        map<pair<unsigned int, unsigned int>, vector<unsigned int>> LatiFacce;
        for (unsigned int i = 0; i < numero_facce; i++) 
        {
            const vector<unsigned int>& faccia = P.Cell2DsVertices[i];
            unsigned int n = faccia.size();
            for (unsigned int j = 0; j < n; j++) 
            {
                unsigned int v1 = faccia[j];
                unsigned int v2 = faccia[(j + 1) % n];
                auto lato = minmax(v1, v2);
                LatiFacce[lato].push_back(i);
            }
        }

        // Per ogni coppia di facce adiacenti genero un lato nel duale
        unsigned int id_lato = 0;
        vector<Vector2i> tmp_lati;

        for (const auto& [lato, facce] : LatiFacce) 
        {
            if (facce.size() == 2) 
            {
                unsigned int f1 = facce[0];
                unsigned int f2 = facce[1];
                P_duale.Cell1DsId.push_back(id_lato++);
                tmp_lati.push_back(Vector2i(f1, f2));
            }
        }

        P_duale.NumCell1Ds = id_lato;
        P_duale.Cell1DsExtrema.resize(2, id_lato);
        for (unsigned int i = 0; i < id_lato; i++)
        {
            P_duale.Cell1DsExtrema.col(i) = tmp_lati[i];
        }


        // Ogni vertice originale diventa una faccia nel duale
        map<unsigned int, vector<unsigned int>> VerticiFacce;
        for (unsigned int i = 0; i < numero_facce; i++)
        {
            for (unsigned int vid : P.Cell2DsVertices[i])
            {
                VerticiFacce[vid].push_back(i);
            }
        }

        unsigned int id_faccia = 0;
        for (const auto& [vid, facce] : VerticiFacce)
        {
            vector<unsigned int> faccia = facce;
            P_duale.Cell2DsId.push_back(id_faccia);
            P_duale.Cell2DsVertices.push_back(faccia);
            ++id_faccia;
        }
        P_duale.NumCell2Ds = id_faccia;

        // Assegno la cella 3D del duale
        P_duale.NumCell3Ds = 1;
        P_duale.Cell3DsFaces.resize(1);
        for (unsigned int i = 0; i < P_duale.NumCell2Ds; i++) 
        {
            P_duale.Cell3DsFaces.push_back(i);
        }

        P_duale.Cell3DsVertices = P_duale.Cell0DsId;
        P_duale.Cell3DsEdges = P_duale.Cell1DsId;
        P_duale.Cell3DsFaces = P_duale.Cell2DsId;

        // Costruisco della mappa dei lati nel duale
        map<pair<unsigned int, unsigned int>, unsigned int> MappaLati;
        for (unsigned int eid = 0; eid < P_duale.NumCell1Ds; eid++)
        {
            unsigned int v1 = P_duale.Cell1DsExtrema(0, eid);
            unsigned int v2 = P_duale.Cell1DsExtrema(1, eid);
            auto key = minmax(v1, v2);
            MappaLati[key] = eid;
        }

        // Riempimento Cell2DsEdges
        P_duale.Cell2DsEdges.resize(P_duale.NumCell2Ds);
        for (unsigned int fid = 0; fid < P_duale.NumCell2Ds; fid++) 
        {
            const vector<unsigned int>& vertici = P_duale.Cell2DsVertices[fid];
            vector<unsigned int> id_lati;
            unsigned int n = vertici.size();

            for (unsigned int i = 0; i < n; i++) 
            {
                unsigned int v1 = vertici[i];
                unsigned int v2 = vertici[(i + 1) % n];
                auto key = minmax(v1, v2);

                auto it = MappaLati.find(key);
                if (it != MappaLati.end()) 
                {
                    id_lati.push_back(it->second);
                } 
                else 
                {   
                    MappaLati[key] = id_lato++;
                }
                id_lati.push_back(MappaLati[key]);
            }
            // Popola Cell2DsEdges per ogni faccia del duale
            P_duale.Cell2DsEdges[fid] = id_lati;
        }

        return P_duale;   
    }


    // ***************************************************************************
    // FUNZIONE ProiezioneSullaSfera
    // Calcolo la proiezione radiale di ogni suo vertice sulla sfera unitaria centrata nell’origine
    Polyhedron ProiezioneSullaSfera(const Polyhedron& P)
    {
        Polyhedron proiettato = P;
        // Scorro tutte le colonne della matrice Cell0DsCoordinates (vertici del poliedro) e li normalizzo
        for (int i = 0; i < proiettato.Cell0DsCoordinates.cols(); i++)
        {
            proiettato.Cell0DsCoordinates.col(i).normalize(); 
        }

        return proiettato;
    }


    // ***************************************************************************
    // FUNZIONE EsportazionePoliedro
    // Esporta i vertici e lati in file .inp per la visualizzazione su Paraview
    void EsportazionePoliedro(const Polyhedron& P, const vector<Gedim::UCDProperty<double>>& points_properties, const vector<Gedim::UCDProperty<double>>& segments_properties)
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
        
        /*
        {
            utilities.ExportPolygons("./Cell2Ds.inp",
                P.Cell0DsCoordinates,
                P.Cell2DsVertices,
                points_properties,
                segments_properties);
        }
        */
       
        EsportazioneFile(P);
    }


    // ***************************************************************************
    // FUNZIONE EsportazioneFile
    // Esporta su file .txt tutte le informazioni di un oggetto Polyhedron   
    void EsportazioneFile(const Polyhedron& P)
    {
        // File: "Cell0Ds.txt"
        ofstream output_file0("Cell0Ds.txt");
        output_file0 << "Id;X;Y;Z\n"; 
        
        const MatrixXd& scrivi_vertici = P.Cell0DsCoordinates;
        for(unsigned int id = 0; id < P.NumCell0Ds; id++)
            output_file0 << defaultfloat << id << ';' << scientific << setprecision(16) << scrivi_vertici(0,id) << ';' << scrivi_vertici(1,id) << ';' << scrivi_vertici(2,id) << '\n';
        output_file0.close();

        // File: "Cell1Ds.txt"
        ofstream output_file1("Cell1Ds.txt");
        output_file1 << "Id;Origin;End\n"; // header

        const MatrixXi& scrivi_lati = P.Cell1DsExtrema;
        for(unsigned int id = 0; id < P.NumCell1Ds; id++)
            output_file1 << id << ';' << scrivi_lati(0,id) << ';' << scrivi_lati(1,id) << '\n';
        output_file1.close();

        // File: "Cell2Ds.txt"
        ofstream output_file2("Cell2Ds.txt");
        output_file2 << "Id;NumVertices;Vertices;NumEdges;Edges\n"; 
        
        const vector<vector<unsigned int>>& scrivi_vertici_facce = P.Cell2DsVertices;
        const vector<vector<unsigned int>>& scrivi_lati_facce = P.Cell2DsEdges;
        for(unsigned int id = 0; id < P.NumCell2Ds; id++)
        {
            output_file2 << id << ';' << scrivi_vertici_facce[id].size();
            for(unsigned int i = 0; i < scrivi_vertici_facce[id].size(); i++)
                output_file2 << ';' << scrivi_vertici_facce[id][i];
            
            output_file2 << ';' << scrivi_lati_facce[id].size();
            for(unsigned int j = 0; j < scrivi_lati_facce[id].size(); j++)
                output_file2 << ';' << scrivi_lati_facce[id][j];

            output_file2 << '\n';
        }
        output_file2.close();

        // File: "Cell3Ds.txt"
        ofstream output_file3("Cell3Ds.txt");
        output_file3 << "Id;NumVertices;Vertices;NumEdges;Edges;NumFaces;Faces\n";

        const vector<unsigned int>& scrivi_vertici_poliedro = P.Cell3DsVertices;
        const vector<unsigned int>& scrivi_lati_poliedro = P.Cell3DsEdges;
        const vector<unsigned int>& scrivi_facce_poliedro = P.Cell3DsFaces;

        // visto che abbiamo un solo poliedro, riscriviamo direttamente i dati della struct senza il ciclo
        unsigned int id =  P.Cell3DsId;
        output_file3 << id << ';' << scrivi_vertici_poliedro.size();
        for(unsigned int i = 0; i < scrivi_vertici_poliedro.size(); i++)
            output_file3 << ';' << scrivi_vertici_poliedro[i];
        
        output_file3 << ';' << scrivi_lati_poliedro.size();
        for(unsigned int j = 0; j < scrivi_lati_poliedro.size(); j++)
            output_file3 << ';' << scrivi_lati_poliedro[j];

        output_file3 << ';' << scrivi_facce_poliedro.size();
        for(unsigned int j = 0; j < scrivi_facce_poliedro.size(); j++)
            output_file3 << ';' << scrivi_facce_poliedro[j];

        output_file3 << '\n';

        output_file3.close();
    }


    // ***************************************************************************
    // FUNZIONE CalcoloCamminoMinimo
    // Calcolo il percorso più breve tra due vertici v1 e v2 all'interno del poliedro P
    pair<vector<unsigned int>, double> CalcoloCamminoMinimo(const Polyhedron& P, const unsigned int& v1, const unsigned int& v2)
    {
        const MatrixXi& lati = P.Cell1DsExtrema;
        const MatrixXd& coordinate = P.Cell0DsCoordinates;

        // Costruisco lista adiacenza con pesi (grafo non orientato)
        vector<vector<pair<unsigned int, double>>> lista_adiacenza(P.NumCell0Ds); 
        for(unsigned int j = 0; j < lati.cols(); j++)
        {   
            // Estraggo gli estremi del lato j
            unsigned int inizio = lati(0, j);
            unsigned int fine = lati(1, j);
            // Calcolo la distanza euclidea tra i due vertici
            Vector3d p1 = coordinate.col(inizio);
            Vector3d p2 = coordinate.col(fine);
            double distanza = (p1 - p2).norm();
            // Aggiungo archi in entrambe le direzioni (grafo non orientato)
            lista_adiacenza[inizio].emplace_back(fine, distanza);
            lista_adiacenza[fine].emplace_back(inizio, distanza);
        }

        // Implemento Algoritmo Dijkstra
        const unsigned int& N = P.NumCell0Ds;
        // Array delle distanze minime inizializzate a infinito
        vector<double> distanza(N, numeric_limits<double>::max());
        // Array dei predecessori per ricostruire il cammino
        vector<unsigned int> parent(N, N); 
        // Coda prioritaria (min-heap) con la distanza minima come chiave
        priority_queue<pair<double, unsigned int>, vector<pair<double, unsigned int>>, greater<pair<double, unsigned int>>> heap;

        // Il nodo iniziale ha distanza zero
        distanza[v1] = 0.0;
        heap.emplace(0.0, v1);

        // Ciclo principale di Dijkstra
        while(!heap.empty())
        {
            auto [dist_corrente, nodo_corrente] = heap.top();
            heap.pop();

            // Se raggiunto nodo finale, termina
            if(nodo_corrente == v2) 
                break; 
            
            // Esplora tutti i vicini del nodo corrente
            for(auto [vicino, lunghezza] : lista_adiacenza[nodo_corrente])
            {
                double nuova_distanza = dist_corrente + lunghezza;
                // Se trovo un cammino più corto verso "vicino", aggiorno
                if(nuova_distanza < distanza[vicino])
                {
                    distanza[vicino] = nuova_distanza;
                    parent[vicino] = nodo_corrente;
                    heap.emplace(nuova_distanza, vicino);
                }
            }
        }

        // Ricostruisco il cammino minimo
        vector<unsigned int> cammino;
        // Se la distanza finale è ancora infinito, non c'è percorso
        if(distanza[v2] == numeric_limits<double>::max())
        {
            cout << "Nessun percorso trovato tra " << v1 << " e " << v2 << endl;
            return {cammino, distanza[v2]};
        }

        // Risalgo il cammino partendo da v2 fino a v1 usando parent[]
        for(unsigned int nodo = v2; nodo != N; nodo = parent[nodo])
        {
            cammino.push_back(nodo);
        }
        // Inverto il cammino perché è stato costruito da fine a inizio
        reverse(cammino.begin(), cammino.end());

        return {cammino, distanza[v2]};
    }

}

