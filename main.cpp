#include <iostream>
#include <vector>
#include "Polyhedron.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"
#include <unordered_set>

using namespace std;
using namespace Eigen;
using namespace PolyhedronLibrary;

int main(int argc, char* argv[])
{
	unsigned int p = 0;
	unsigned int q = 0;
	unsigned int b = 0;
	unsigned int c = 0;
	unsigned int v1 = 0;
	unsigned int v2 = 0;
	
	// controllo che il tipo dei valori inseriti sia int
	p = stoi(argv[1]); //commenta stoi
	q = stoi(argv[2]);
	b = stoi(argv[3]);
	c = stoi(argv[4]);
	if (argc == 5) 
	{
		cout << "Input ricevuto: p = " << p << ", q = " << q << ", b = " << b << ", c = " << c << "\n";
    }
	else if (argc == 7) 
	{
		v1 = stoi(argv[5]);
		v2 = stoi(argv[6]);
		cout << "Input ricevuto: p = " << p << ", q = " << q << ", b = " << b << ", c = " << c << ", v1 = " << v1 << ", v2 = " << v2 << "\n";
    }
	else
 	{
		cerr << "Inserimento non valido. Inserire 4 o 6 numeri interi positivi.\n" <<endl;
		return 1;
	}

	// controllo le condizioni su p e q
	if (p < 3 || q < 3 || double (1.0/p)+double (1.0/q) <= 0.5)
	{
		cerr << "Errore: i parametri p e q non rispettano le condizioni dei poligoni platonici.\n";
    	return 1;
	}

	// controllo delle classi di triangolazione valide
	if (b != c && b!= 0 && c!= 0)
	{
		cerr << "Errore: classe non trattata.\n";
		return 1;
	}

	// Creo il solido platonico e riempo la struct
	// CASO TETRAEDRO, OTTAEDRO, ICOSAEDRO
	Polyhedron Poliedro_finale;

	if (p == 3)
	{	
		Polyhedron mesh = buildPlatonicSolid(q); //definiamo il solido platonico di partenza

		if (b != c) // TRIANGOLAZIONE DI CLASSE I
		{
			cout<<"Eseguiamo la triangolazione di classe I"<<endl;
			unsigned int t_value = b + c; //valore che mi indica in quante parti dividere ogni lato del triangolo
			Polyhedron mesh_triangolata1 = triangulateClass1(mesh, t_value); //triangolazione di classe 1 con parametro t_value
			Polyhedron mesh_proiettata = projectPolyhedronOnSphere(mesh_triangolata1);
			Poliedro_finale = mesh_proiettata;
		}

		else // TRIANGOLAZIONE DI CLASSE II
		{
			cout<<"Eseguiamo la triangolazione di classe II"<<endl;
			Polyhedron mesh_triangolata2 = triangulateClass2(mesh, b);
			Polyhedron mesh_proiettata = projectPolyhedronOnSphere(mesh_triangolata2);
			Poliedro_finale = mesh_proiettata;
		}
	}
	

	// CASO CUBO (DUALE DELL'OTTAEDRO) e CASO DODECAEDRO (DUALE DELL'ICOSAEDRO)
	else if (p == 4 || p == 5)  
	{
		Polyhedron mesh = buildPlatonicSolid(p); //definiamo il solido platonico di partenza

		if (b != c) // TRIANGOLAZIONE DI CLASSE I
		{
			cout<<"Eseguiamo la triangolazione di classe I"<<endl;
			unsigned int t_value = b + c; //valore che mi indica in quante parti dividere ogni lato del triangolo
			Polyhedron mesh_triangolata1 = triangulateClass1(mesh, t_value); //triangolazione di classe 1 con parametro t_value	
			Polyhedron mesh_dualizzata = Dualize(mesh_triangolata1);
			Polyhedron mesh_proiettata = projectPolyhedronOnSphere(mesh_dualizzata);
			Poliedro_finale = mesh_proiettata;
		}
		else // TRIANGOLAZIONE DI CLASSE II
		{
			cout<<"Eseguiamo la triangolazione di classe II"<<endl;
			Polyhedron mesh_triangolata2 = triangulateClass2(mesh, b);
			Polyhedron mesh_dualizzata = Dualize(mesh_triangolata2);
			Polyhedron mesh_proiettata = projectPolyhedronOnSphere(mesh_dualizzata);
			Poliedro_finale = mesh_proiettata;
		}
	}

	
	//Cammini minimi	
	if(argc == 7)
	{	
		// controllo su v1 e v2
		if (v1 >= Poliedro_finale.NumCell0Ds || v2 >= Poliedro_finale.NumCell0Ds) 
		{
			cerr<<"I valori di v1 e v2 non sono validi: devono appartenere all'intervallo [0, " << Poliedro_finale.NumCell0Ds -1<< "]" <<endl;
			return 1;
		}

		vector<Gedim::UCDProperty<double>> points_properties;
		vector<Gedim::UCDProperty<double>> segments_properties;
		points_properties.reserve(1); //Proprietà: ShortPath
		segments_properties.reserve(1);

		vector<double> prop_vert(Poliedro_finale.NumCell0Ds, 0.0);
		vector<double> prop_edges(Poliedro_finale.NumCell1Ds, 0.0);

		//Riempimento point_properties
		Gedim::UCDProperty<double> pointP;
		pointP.Label = "ShortPath";
		pointP.NumComponents = 1;
		pointP.Data = prop_vert.data();
		points_properties.push_back(pointP);

		//Riempimento segments_properties
		Gedim::UCDProperty<double> edgeP;
		edgeP.Label = "ShortPath";
		edgeP.NumComponents = 1;
		edgeP.Data = prop_edges.data();
		segments_properties.push_back(edgeP);


		vector<unsigned int> path = Cammini_minimi(Poliedro_finale, v1, v2);
	
		//Riempimento prop_vert
		for(unsigned int id_v : path) prop_vert[id_v] = 1.0; //imposta a 1.0 la proprietà dei vertici lungo il percorso

		// Riempimento prop_edges
		for(unsigned int i=0; i < path.size() - 1; i++)
		{
			unsigned int corrente = path[i];
			unsigned int successivo = path[i+1];
			std::unordered_set<unsigned int> edge_ref = {corrente, successivo};

			// Itera sugli spigoli e, se lo spigolo è formato da corrente, successivo, imposta prop_edges[j] = 1.0
			for(unsigned int j=0; j < Poliedro_finale.Cell1DsExtrema.cols(); j++)
			{
				Eigen::Vector2i edge = Poliedro_finale.Cell1DsExtrema.col(j);
				unsigned int a = edge[0]; 
				unsigned int b = edge[1];
				if(edge_ref.count(a) + edge_ref.count(b) == 2) prop_edges[j] = 1.0;
			}
		}

		ExportPolyhedron(Poliedro_finale, points_properties, segments_properties);
	}

	else
	{
		ExportPolyhedron(Poliedro_finale);
	}
	cout << "I dati del solido geodetico sono stati salvati nei file: 'Cell0Ds.txt','Cell1Ds.txt','Cell2Ds.txt','Cell3Ds.txt' " << endl;
	
	return 0;

}


