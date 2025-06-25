#pragma once  // assicura che i file vengano inclusi solo una volta per unità di compilazione

#include "UCDUtilities.hpp"
#include "Polyhedron.hpp"

#include <iostream>

using namespace std;

namespace PolyhedronLibrary
{
    // La funzione SalvataggioDati serve per costruire i dati geometrici (vertici e facce) di un solido platonico, 
    // dati in input un intero q che rappresenta il tipo di solido (q= 3 tetraedro, q=4 ottaedro, q=5 icosaedro).
    // Restituisce una coppia contenente: le coordinate normalizzate dei vertici del solido, le facce espresse come 
    // insiemi di id dei vertici. 
    pair<vector<Vector3d>, vector<vector<unsigned int>>> SalvataggioDati(unsigned int& q);


    // La funzione CreazioneSolidoPlatonico costruisce e restituisce un oggetto Polyhedron, che rappresenta un solido 
    // platonico tridimensionale specificato da q. Organizza i dati in celle di dimensione 0 (vertici), 1 (lati), 
    // 2 (facce) e 3 (il solido stesso).
    Polyhedron CreazioneSolidoPlatonico(unsigned int& q);


    // La funzione TriangolazioneClasse1 prende un poliedro P con facce triangolari e restituisce un nuovo poliedro P1 
    // in cui ogni faccia è ulteriormente suddivisa in triangoli più piccoli tramite una triangolazione detta “di classe 1”, 
    // regolata da un parametro intero t_value. Restituisce un elemento Polyhedron con i dati della triangolazione 1.
    Polyhedron TriangolazioneClasse1(const Polyhedron& P, const unsigned int& t_value);


    // La funzione TriangolazioneClasse2 effettua una raffinazione di classe 2 su un poliedro, basandosi sulla triangolazione 
    // di classe 1 e sul valore in input b: aggiunge baricentri delle facce e punti medi dei lati, collegandoli per formare 
    // nuovi triangoli. Restituisce un elemento Polyhedron con i dati della triangolazione 2.
    Polyhedron TriangolazioneClasse2(const Polyhedron& P, const unsigned int& b);


    // La funzione Dualizzazione costruisce il duale di un poliedro dato in input, cioè un nuovo poliedro dove: ogni faccia del
    // poliedro originale diventa un vertice del duale, ogni vertice del poliedro originale diventa una faccia del duale, ogni 
    // lato tra due facce adiacenti nel poliedro originale diventa un lato nel duale che collega i relativi nuovi vertici 
    // (baricentri delle facce). Restituisce un elemento Polyhedron con i dati della dualizzazione.
    Polyhedron Dualizzazione(const Polyhedron& P);


    // La funzione ProiezioneSullaSfera proietta i vertici di un poliedro sulla superficie della sfera unitaria, cioè 
    // normalizzati in modo che abbiano tutti distanza 1 dall’origine. Restituisce un elemento Polyhedron con i punti normalizzati.
    Polyhedron ProiezioneSullaSfera(const Polyhedron& P);


    // La funzione EsportazionePoliedro esporta i punti e lati di Polyhedron in file esterni .inp, usando la libreria 
    // Gedim::UCDUtilities. È una funzione void (non restituisce nulla).
    void EsportazionePoliedro(const Polyhedron& P,
                            const vector<Gedim::UCDProperty<double>>& points_properties = {},
                            const vector<Gedim::UCDProperty<double>>& segments_properties = {});


    // La funzione EsportazioneFile esporta su file .txt tutte le informazioni strutturali di un oggetto Polyhedron:
    // i vertici (0D) → Cell0Ds.txt, i lati (1D) → Cell1Ds.txt, le facce (2D) → Cell2Ds.txt, le celle (3D) → Cell3Ds.txt.
    // È una funzione void (non restituisce nulla).
    void EsportazioneFile(const Polyhedron& P);


    // La funzione CalcoloCamminoMinimo calcola il percorso più breve tra due vertici v1 e v2 all'interno del poliedro P, 
    // usando l'algoritmo di Dijkstra, dove il peso degli archi è la distanza euclidea tra vertici connessi.
    // Restituisce un vector con gli ID dei vertici nel cammino minimo.
    pair<vector<unsigned int>, double> CalcoloCamminoMinimo(const Polyhedron& P, const unsigned int& v1, const unsigned int& v2);
}
