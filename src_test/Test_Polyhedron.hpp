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
        // Crea il numero minimo di elementi salvabili sui file .inp (2 punti e 1 segmento)
        Polyhedron polyhedron;

        // Aggiungi due vertici
        polyhedron.Cell0DsCoordinates.resize(3, 2);
        polyhedron.Cell0DsCoordinates.col(0) << 0.0, 0.0, 0.0;
        polyhedron.Cell0DsCoordinates.col(1) << 1.0, 0.0, 0.0;
        polyhedron.Cell0DsId = {0, 1};
        polyhedron.NumCell0Ds = 2;

        // Aggiungi un segmento che collega i due vertici
        polyhedron.Cell1DsExtrema.resize(2, 1);
        polyhedron.Cell1DsExtrema(0, 0) = 0;
        polyhedron.Cell1DsExtrema(1, 0) = 1;
        polyhedron.Cell1DsId = {0};
        polyhedron.NumCell1Ds = 1;

        polyhedron.Cell2DsId = {};
        polyhedron.NumCell2Ds = 0;
        polyhedron.Cell3DsId = {};
        polyhedron.NumCell3Ds = 0;

        // Chiama ExportPolyhedron SENZA proprietà
        ExportPolyhedron(polyhedron, {}, {});

        // Verifica che i file .inp siano stati creati e non siano vuoti
        std::ifstream file0("Cell0Ds.inp");
        ASSERT_TRUE(file0.is_open()) << "Cell0Ds.inp non è stato creato";
        ASSERT_NE(file0.peek(), std::ifstream::traits_type::eof()) << "Cell0Ds.inp è vuoto";
        // .peek() = prende il primo carattere del file, NE = not equal
        file0.close();

        std::ifstream file1("Cell1Ds.inp");
        ASSERT_TRUE(file1.is_open()) << "Cell1Ds.inp non è stato creato";
        ASSERT_NE(file1.peek(), std::ifstream::traits_type::eof()) << "Cell1Ds.inp è vuoto";
        file1.close();

        // Svota i file generati
        std::remove("Cell0Ds.inp");
        std::remove("Cell1Ds.inp");
    }





}