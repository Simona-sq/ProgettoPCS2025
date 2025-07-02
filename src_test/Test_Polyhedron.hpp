#pragma once

#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include <fstream>
#include "Utils.hpp"
#include "Eigen/Eigen"
#include "Polyhedron.hpp"
#include "Utils.hpp"


using namespace Eigen;
using namespace PolyhedronLibrary;

namespace PolyhedronLibrary {

    // G TEST SU SalvataggioDati
    TEST(TestPolyhedron, TestSalvataggioDati)
    {
        unsigned int q = 5;  // Test per l'icosaedro
        auto [vertici, facce] = SalvataggioDati(q);

        for (const auto& faccia : facce)
            EXPECT_EQ(faccia.size(), 3);

        for (const auto& vertice : vertici){
            EXPECT_NEAR(vertice.norm(), 1.0, 2.2e-16);}
    }


    // G TEST SU CreazioneSolidoPlatonico
    TEST(TestPolyhedron, TestCreazioneSolidoPlatonico)
    {
        unsigned int q = 5;  // Test per l'icosaedro
        Polyhedron P = CreazioneSolidoPlatonico(q);

        EXPECT_EQ(P.NumCell0Ds, 12);
        EXPECT_EQ(P.NumCell1Ds, 30);
        EXPECT_EQ(P.NumCell2Ds, 20);
        EXPECT_EQ(P.NumCell3Ds, 1);
        EXPECT_EQ(P.Cell0DsCoordinates.cols(), 12);
    }


    // G TEST su TriangolazioneClasse1
    TEST(TestPolyhedron, TestTriangolazioneClasse1)
    {
        unsigned int q = 3;  // Test per il tetraedro
        unsigned int t = 2; 
        unsigned int T = t * t;  // T = b^2 + bc + c^2 con t = b + c ⇒ T = t^2

        // Dati delle formule
        unsigned int V_atteso = 2 * T + 2;
        unsigned int E_atteso = 6 * T;
        unsigned int F_atteso = 4 * T;

        // Costruzione solido e triangolazione
        Polyhedron P = CreazioneSolidoPlatonico(q);
        Polyhedron P1 = TriangolazioneClasse1(P, t);

        // Verifica numero vertici
        EXPECT_EQ(P1.NumCell0Ds, V_atteso) << "Numero vertici errato";

        // Verifica numero lati
        EXPECT_EQ(P1.NumCell1Ds, E_atteso) << "Numero lati errato";

        // Verifica numero facce
        EXPECT_EQ(P1.NumCell2Ds, F_atteso) << "Numero facce errato";
    }


    // G TEST su TriangolazioneClasse2
    TEST(TestPolyhedron, TestTriangolazioneClasse2)
    {
        unsigned int q = 3;  // Test per il tetraedro
        unsigned int b = 2;

        // Costruisco il solido platonico 
        Polyhedron P = CreazioneSolidoPlatonico(q);

        // Triangolazione di classe 2
        Polyhedron P2 = TriangolazioneClasse2(P, b);

        // Estrai il numero originale di vertici, lati e facce
        unsigned int V0 = P.NumCell0Ds; // 4
        unsigned int E0 = P.NumCell1Ds; // 6
        unsigned int F0 = P.NumCell2Ds; // 4

        // Calcolo delle quantità attese secondo le formule della lavagna
        unsigned int V_atteso = V0  + E0 * (2 * b - 1) + F0 * ((3 * b * b) / 2 - (3 * b) / 2 + 1);

        unsigned int E_atteso = E0 * (2 * b) + F0 * ((9 * b * b) / 2 + (3 * b) / 2);

        unsigned int F_atteso = F0 * (3 * b * b + 3 * b);

        EXPECT_EQ(P2.NumCell0Ds, V_atteso) << "Numero di vertici errato";
        EXPECT_EQ(P2.NumCell1Ds, E_atteso) << "Numero di lati errato";
        EXPECT_EQ(P2.NumCell2Ds, F_atteso) << "Numero di facce errato";
    }


    // G TEST SU Dualizzazione
    TEST(TestPolyhedron, TestDualizzazione)
    {
        unsigned int q = 3;   // Test per il tetraedro
        Polyhedron P = CreazioneSolidoPlatonico(q);

        // Dualizza il tetraedro
        Polyhedron P_duale = Dualizzazione(P);

        // Il duale di un tetraedro ha:
        // - 4 facce nel solido di partenza → 4 vertici nel duale
        // - 6 spigoli nel solido di partenza → 6 spigoli nel duale
        // - 4 vertici nel solido di partenza → 4 facce nel duale

        unsigned int V_atteso = P.NumCell2Ds;
        unsigned int E_atteso = P.NumCell1Ds;
        unsigned int F_atteso = P.NumCell0Ds;

        EXPECT_EQ(P_duale.NumCell0Ds, V_atteso) << "Numero di vertici errato";
        EXPECT_EQ(P_duale.NumCell1Ds, E_atteso) << "Numero di lati errato";
        EXPECT_EQ(P_duale.NumCell2Ds, F_atteso) << "Numero di facce errato";
    }


    // G TEST SU ProiezioneSullaSfera
    TEST(TestPolyhedron, TestProiezioneSullaSfera)
    {
        Polyhedron P;
        // Inserisco 3 vertici arbitrari NON normalizzati
        P.NumCell0Ds = 3;
        P.Cell0DsId = {0, 1, 2};
        P.Cell0DsCoordinates.resize(3, 3);
        P.Cell0DsCoordinates.col(0) << 1.0, 2.0, 2.0; 
        P.Cell0DsCoordinates.col(1) << 0.0, 3.0, 4.0;
        P.Cell0DsCoordinates.col(2) << -1.0, -1.0, -1.0;

        // Eseguo proiezione con la nostra funzione
        Polyhedron P_proiettato = ProiezioneSullaSfera(P);

        // Verifico che la dimensione resti la stessa
        ASSERT_EQ(P_proiettato.NumCell0Ds, 3);
        ASSERT_EQ(P_proiettato.Cell0DsCoordinates.cols(), 3);

        // Verifico che ogni punto abbia norma ≈ 1
        for (int i = 0; i < P_proiettato.Cell0DsCoordinates.cols(); ++i) 
        {
            double norma = P_proiettato.Cell0DsCoordinates.col(i).norm();
            EXPECT_NEAR(norma, 1.0, 2.2e-16) << "Vertice " << i << " non è unitario, ha norma = " << norma << endl;
        }

        // Verifico che gli ID siano invariati
        EXPECT_EQ(P_proiettato.Cell0DsId, P.Cell0DsId);
    }

    
    // G TEST PER EsportazionePoliedro su un segmento di due punti
    TEST(TestPolyhedron, TestEsportazionePoliedro)
    {
        Polyhedron P;

        // Vertici
        P.Cell0DsCoordinates.resize(3, 2);
        P.Cell0DsCoordinates.col(0) << 0.0, 0.0, 0.0;
        P.Cell0DsCoordinates.col(1) << 1.0, 0.0, 0.0;
        P.Cell0DsId = {0, 1};
        P.NumCell0Ds = 2;

        // Segmento
        P.Cell1DsExtrema.resize(2, 1);
        P.Cell1DsExtrema(0, 0) = 0;
        P.Cell1DsExtrema(1, 0) = 1;
        P.Cell1DsId = {0};
        P.NumCell1Ds = 1;

        // Facce (vuote)
        P.Cell2DsVertices = {};
        P.Cell2DsEdges = {};
        P.Cell2DsId = {};
        P.NumCell2Ds = 0;

        // Poliedro 3D 
        P.Cell3DsVertices = {{0, 1}};
        P.Cell3DsEdges = {{0}};
        P.Cell3DsFaces = {{}};
        P.Cell3DsId = 0;
        P.NumCell3Ds = 1;

        // Esportazione
        EsportazionePoliedro(P, {}, {});

        // Controlli su file
        ifstream file0("Cell0Ds.inp");
        ASSERT_TRUE(file0.is_open()) << "Cell0Ds.inp non è stato creato";

        ASSERT_NE(file0.peek(), ifstream::traits_type::eof()) << "Cell0Ds.inp è vuoto";  // .peek() = prende il primo carattere del file, NE = not equal
        file0.close();

        ifstream file1("Cell1Ds.inp");
        ASSERT_TRUE(file1.is_open()) << "File non aperto";
        ASSERT_NE(file1.peek(), -1) << "File vuoto";
        file1.close();

        // Pulizia
        remove("Cell0Ds.inp");
        remove("Cell1Ds.inp");
        remove("Cell0Ds.txt");
        remove("Cell1Ds.txt");
        remove("Cell2Ds.txt");
        remove("Cell3Ds.txt");
    }
    

    // G TEST PER Esporta_file su un segmento di due punti
    TEST(TestPolyhedron, TestEsporta_file) 
    {
        // Creo un polyhedron minimale: 2 vertici, 1 lato, 1 faccia, 1 cella
        Polyhedron P;
    
        // Vertici
        P.Cell0DsCoordinates.resize(3, 2);
        P.Cell0DsCoordinates.col(0) << 0.0, 0.0, 0.0;
        P.Cell0DsCoordinates.col(1) << 1.0, 0.0, 0.0;
        P.Cell0DsId = {0, 1};
        P.NumCell0Ds = 2;
    
        // Lato
        P.Cell1DsExtrema.resize(2, 1);
        P.Cell1DsExtrema(0, 0) = 0;
        P.Cell1DsExtrema(1, 0) = 1;
        P.Cell1DsId = {0};
        P.NumCell1Ds = 1;
    
        // Faccia
        P.Cell2DsVertices = {{0, 1}};
        P.Cell2DsEdges = {{0}};
        P.Cell2DsId = {0};
        P.NumCell2Ds = 1;
    
        // Cella 3D
        P.Cell3DsVertices = {{0, 1}};
        P.Cell3DsEdges = {{0}};
        P.Cell3DsFaces = {{0}};
        P.Cell3DsId = 0;
        P.NumCell3Ds = 1;
    
        EsportazioneFile(P);
    
        // Verifico esistenza e non-vuotezza dei file
        const vector<string> files = { "Cell0Ds.txt", "Cell1Ds.txt", "Cell2Ds.txt", "Cell3Ds.txt" };
    
        for (const auto& nome_file : files)
        {
            ifstream file(nome_file);
            ASSERT_TRUE(file.is_open()) << "Il file " << nome_file << " non è stato creato";
            ASSERT_NE(file.peek(), -1) << "File vuoto";
            file.close();
        }
    
        // Pulizia
        for (const auto& nome_file : files)
        {
            remove(nome_file.c_str());
        }
    }

    // G TEST PER CalcoloCamminoMinimo
    TEST(TestPolyhedron, TestCalcoloCamminoMinimo)
    {
        // Testo il cammino minimo tra i vertici 10 e 58 del dodecaedro con triangolazione di classe 1 con t=4
        unsigned int q = 4;
        unsigned int t = 4;
        Polyhedron P = CreazioneSolidoPlatonico(q);
        Polyhedron P1 = TriangolazioneClasse1(P, t);
        Polyhedron P_duale = Dualizzazione(P1);
        Polyhedron P_proiettato = ProiezioneSullaSfera(P_duale);
    
        unsigned int v1 = 10;
        unsigned int v2 = 58;
    
        // Trovo il cammino minimo
        auto [cammino, distanza] = CalcoloCamminoMinimo(P_proiettato, v1, v2);
    
        // Verifiche sul cammino
        ASSERT_FALSE(cammino.empty()) << "Cammino non trovato.";
        EXPECT_EQ(cammino.size() - 1, 3) << "Il cammino deve attraversare 3 lati (4 vertici).";

        // Verifico che la sequenza di vertici sia quella attesa (ipotizziamo un cammino noto trovato con Paraview)
        vector<unsigned int> cammino_atteso = {10, 9, 57, 58}; 
        EXPECT_EQ(cammino, cammino_atteso) << "La sequenza del cammino minimo non è quella attesa";
    }
}