#pragma once

#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include <fstream>
//#include <cstdio>

#include "Utils.hpp"
#include "Eigen/Eigen"
#include "Polyhedron.hpp"
#include "Utils.hpp"


using namespace Eigen;
using namespace PolyhedronLibrary;

namespace PolyhedronLibrary {

    // G TEST SU GetSolidData
    TEST(GetSolidDataTest, ReturnsCorrectStructureForIcosahedron) {
        unsigned int q = 5;  // Test per l'icosaedro
        auto [verts, faces] = getSolidData(q);

        EXPECT_EQ(verts.size(), 12);
        EXPECT_EQ(faces.size(), 20);

        for (const auto& face : faces)
            EXPECT_EQ(face.size(), 3);

        for (const auto& v : verts)
            EXPECT_NEAR(v.norm(), 1.0, 1e-8);
    }


    // G TEST SU buildPlatonicSolid
    TEST(BuildPlatonicSolidTest, IcosahedronStructure) {
        unsigned int q = 5;  // Test per l'icosaedro
        Polyhedron P = buildPlatonicSolid(q);

        EXPECT_EQ(P.NumCell0Ds, 12);
        EXPECT_EQ(P.NumCell1Ds, 30);
        EXPECT_EQ(P.NumCell2Ds, 20);
        EXPECT_EQ(P.NumCell3Ds, 1);
        EXPECT_EQ(P.Cell0DsCoordinates.cols(), 12);

        for (int i = 0; i < P.Cell0DsCoordinates.cols(); ++i) {
            double norm = P.Cell0DsCoordinates.col(i).norm();
            EXPECT_NEAR(norm, 1.0, 1e-8);
        }
    }


    // G TEST su triangulateClass1
    TEST(TestTriangulateClass1, CountsAndValenceTetrahedronT2)
    {
        unsigned int q = 3;  // Test per il tetraedro
        unsigned int t = 2; 
        unsigned int T = t * t;     // T = b^2 + bc + c^2 con t = b + c ⇒ T = t^2

        // Dati delle formule
        unsigned int expected_V = 2 * T + 2;
        unsigned int expected_E = 6 * T;
        unsigned int expected_F = 4 * T;

        // Costruzione solido e triangolazione
        Polyhedron P = buildPlatonicSolid(q);
        Polyhedron P1 = triangulateClass1(P, t);

        // Verifica numero vertici
        EXPECT_EQ(P1.NumCell0Ds, expected_V) << "Numero vertici errato";

        // Verifica numero lati
        EXPECT_EQ(P1.NumCell1Ds, expected_E) << "Numero lati errato";

        // Verifica numero facce
        EXPECT_EQ(P1.NumCell2Ds, expected_F) << "Numero facce errato";
    }


    // G TEST su triangulateClass2
    TEST(TestTriangulateClass2, TetrahedronB2Counts)
    {
        unsigned int q = 3;  // Test per il tetraedro
        unsigned int b = 2;

        // Costruisci il solido platonico (tetraedro)
        Polyhedron P = buildPlatonicSolid(q);

        // Triangolazione di classe 2
        Polyhedron P2 = triangulateClass2(P, b);

        // Estrai il numero originale di vertici, lati e facce
        unsigned int V0 = P.NumCell0Ds; // 4
        unsigned int E0 = P.NumCell1Ds; // 6
        unsigned int F0 = P.NumCell2Ds; // 4

        // Calcolo delle quantità attese secondo le formule della lavagna
        unsigned int expectedV = V0  + E0 * (2 * b - 1) + F0 * ((3 * b * b) / 2 - (3 * b) / 2 + 1);

        unsigned int expectedE = E0 * (2 * b) + F0 * ((9 * b * b) / 2 + (3 * b) / 2);

        unsigned int expectedF = F0 * (3 * b * b + 3 * b);

        // Asserzioni
        EXPECT_EQ(P2.NumCell0Ds, expectedV) << "Numero di vertici errato";
        EXPECT_EQ(P2.NumCell1Ds, expectedE) << "Numero di lati errato";
        EXPECT_EQ(P2.NumCell2Ds, expectedF) << "Numero di facce errato";
    }


    // G TEST SU Dualize
    TEST(TestDualize, CheckDualOnTetrahedron)
    {
        unsigned int q = 3;   // Test per il tetraedro
        Polyhedron P = buildPlatonicSolid(q);

        // Dualizza il tetraedro
        Polyhedron P_duale = Dualize(P);

        // Il duale di un tetraedro ha:
        // - 4 facce nel solido di partenza → 4 vertici nel duale
        // - 6 spigoli nel solido di partenza → 6 spigoli nel duale
        // - 4 vertici nel solido di partenza → 4 facce nel duale
        unsigned int expectedV = P.NumCell2Ds;
        unsigned int expectedE = P.NumCell1Ds;
        unsigned int expectedF = P.NumCell0Ds;

        EXPECT_EQ(P_duale.NumCell0Ds, expectedV) << "Numero di vertici errato";
        EXPECT_EQ(P_duale.NumCell1Ds, expectedE) << "Numero di lati errato";
        EXPECT_EQ(P_duale.NumCell2Ds, expectedF) << "Numero di facce errato";

    }


    // G TEST SU ProjectPolyhedronOnSphere
    TEST(ProjectPolyhedronOnSphereTest, ProjectsVerticesCorrectly)
    {
        Polyhedron P;
        // Inserisci 3 vertici arbitrari NON normalizzati
        P.NumCell0Ds = 3;
        P.Cell0DsId = {0, 1, 2};
        P.Cell0DsCoordinates.resize(3, 3);
        P.Cell0DsCoordinates.col(0) << 1.0, 2.0, 2.0; 
        P.Cell0DsCoordinates.col(1) << 0.0, 3.0, 4.0;
        P.Cell0DsCoordinates.col(2) << -1.0, -1.0, -1.0;

        // Esegui proiezione con la nostra funzione
        Polyhedron projected = projectPolyhedronOnSphere(P);

        // Verifica che la dimensione resti la stessa
        ASSERT_EQ(projected.NumCell0Ds, 3);
        ASSERT_EQ(projected.Cell0DsCoordinates.cols(), 3);

        // Verifica che ogni punto abbia norma ≈ 1
        for (int i = 0; i < projected.Cell0DsCoordinates.cols(); ++i) 
        {
            double norm = projected.Cell0DsCoordinates.col(i).norm();
            EXPECT_NEAR(norm, 1.0, 1e-8) << "Vertice " << i << " non è unitario, ha norma = " << norm << endl;
        }

        // Verifica che gli ID siano invariati
        EXPECT_EQ(projected.Cell0DsId, P.Cell0DsId);
    }


    
    // G TEST PER ExportPolyhedron su un segmento di due punti
    TEST(ExportPolyhedronTest, MinimalSegmentExportWithoutProperties) {
        Polyhedron polyhedron;

        // Vertici
        polyhedron.Cell0DsCoordinates.resize(3, 2);
        polyhedron.Cell0DsCoordinates.col(0) << 0.0, 0.0, 0.0;
        polyhedron.Cell0DsCoordinates.col(1) << 1.0, 0.0, 0.0;
        polyhedron.Cell0DsId = {0, 1};
        polyhedron.NumCell0Ds = 2;

        // Segmento
        polyhedron.Cell1DsExtrema.resize(2, 1);
        polyhedron.Cell1DsExtrema(0, 0) = 0;
        polyhedron.Cell1DsExtrema(1, 0) = 1;
        polyhedron.Cell1DsId = {0};
        polyhedron.NumCell1Ds = 1;

        // Facce (vuote)
        polyhedron.Cell2DsVertices = {};
        polyhedron.Cell2DsEdges = {};
        polyhedron.Cell2DsId = {};
        polyhedron.NumCell2Ds = 0;

        // Poliedro 3D (deve essere inizializzato anche se non ha significato)
        polyhedron.Cell3DsVertices = {{0, 1}};
        polyhedron.Cell3DsEdges = {{0}};
        polyhedron.Cell3DsFaces = {{}};
        polyhedron.Cell3DsId = 0;
        polyhedron.NumCell3Ds = 1;

        // Esportazione
        ExportPolyhedron(polyhedron, {}, {});

        // Controlli su file
        std::ifstream file0("Cell0Ds.inp");
        ASSERT_TRUE(file0.is_open()) << "Cell0Ds.inp non è stato creato";

        // .peek() = prende il primo carattere del file, NE = not equal
        ASSERT_NE(file0.peek(), std::ifstream::traits_type::eof()) << "Cell0Ds.inp è vuoto";
        file0.close();

        std::ifstream file1("Cell1Ds.inp");
        ASSERT_TRUE(file1.is_open()) << "File non aperto";
        ASSERT_NE(file1.peek(), -1) << "File vuoto";
        file1.close();

        // Pulizia
        std::remove("Cell0Ds.inp");
        std::remove("Cell1Ds.inp");
        std::remove("Cell0Ds.txt");
        std::remove("Cell1Ds.txt");
        std::remove("Cell2Ds.txt");
        std::remove("Cell3Ds.txt");
    }
    

    // G TEST PER Esporta_file su un segmento di due punti
    TEST(EsportaFileTest, GeneratesCorrectOutputFiles) {
        // Crea un polyhedron minimale: 2 vertici, 1 lato, 1 faccia, 1 cella
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
    
        // Chiamata alla funzione da testare
        Esporta_file(P);
    
        // Verifica esistenza e non-vuotezza dei file
        const std::vector<std::string> files = {
            "Cell0Ds.txt", "Cell1Ds.txt", "Cell2Ds.txt", "Cell3Ds.txt"
        };
    
        for (const auto& filename : files) {
            std::ifstream file(filename);
            ASSERT_TRUE(file.is_open()) << "Il file " << filename << " non è stato creato";
            ASSERT_NE(file.peek(), -1) << "File vuoto";
            file.close();
        }
    
        // Pulizia finale
        for (const auto& filename : files) {
            std::remove(filename.c_str());
        }
    }

    // G TEST PER cammini_minimi dal vertice 0 al vertice 2
    TEST(CamminiMinimiTest, PathLength5OnTriangulatedTetrahedronClass1T2)
    {
        // 1. Costruzione tetraedro e triangolazione t=2
        unsigned int q = 4;
        Polyhedron P = buildPlatonicSolid(q);
        unsigned int t = 4;
        Polyhedron triangulated = triangulateClass1(P, t);
        Polyhedron dualized = Dualize(triangulated);

        cout << dualized.NumCell0Ds << endl;
    
        unsigned int v_start = 10;
        unsigned int v_end = 58;
    
        // 3. Trova il cammino minimo
        std::vector<unsigned int> path = Cammini_minimi(dualized, v_start, v_end);

        std::cout << "Cammino: ";
        for (auto id : path) std::cout << id << " ";
        std::cout << std::endl;
    
        // 4. Verifiche sul cammino
        ASSERT_FALSE(path.empty()) << "Cammino non trovato.";
        EXPECT_EQ(path.size() - 1, 3) << "Il cammino deve attraversare 3 lati (4 vertici).";

        // 7. Verifica che la sequenza di vertici sia quella attesa (ipotizziamo un cammino noto)
        std::vector<unsigned int> expected_path = {10, 9, 57, 58}; 
        EXPECT_EQ(path, expected_path) << "La sequenza del cammino minimo non è quella attesa";
    }





}