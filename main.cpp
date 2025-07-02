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
	
	p = stoi(argv[1]); // stoi: converte una stringa in un intero
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

	// Controllo le condizioni su p e q
	if (p < 3 || q < 3 || double (1.0/p)+double (1.0/q) <= 0.5)
	{
		cerr << "Errore: i parametri p e q non rispettano le condizioni dei solidi platonici.\n";
    	return 1;
	}

	// Controllo delle classi di triangolazione valide
	if (b != c && b!= 0 && c!= 0)
	{
		cerr << "Errore: classe non trattata.\n";
		return 1;
	}

	// Creo il solido platonico e riempo la struct
	Polyhedron poliedro_finale;

	// CASO TETRAEDRO, OTTAEDRO, ICOSAEDRO
	if (p == 3)
	{	
		Polyhedron mesh = CreazioneSolidoPlatonico(q); // definisco il solido platonico di partenza

		if (b != c) // TRIANGOLAZIONE DI CLASSE I
		{
			cout<<"Eseguiamo la triangolazione di classe I"<<endl;
			unsigned int t_value = b + c; // valore che mi indica in quante parti dividere ogni lato del triangolo
			Polyhedron mesh_triangolata1 = TriangolazioneClasse1(mesh, t_value);
			Polyhedron mesh_proiettata = ProiezioneSullaSfera(mesh_triangolata1);
			poliedro_finale = mesh_proiettata;
		}

		else // TRIANGOLAZIONE DI CLASSE II
		{
			cout<<"Eseguiamo la triangolazione di classe II"<<endl;
			Polyhedron mesh_triangolata2 = TriangolazioneClasse2(mesh, b);
			Polyhedron mesh_proiettata = ProiezioneSullaSfera(mesh_triangolata2);
			poliedro_finale = mesh_proiettata;
		}
	}
	
	// CASO CUBO (DUALE DELL'OTTAEDRO) e CASO DODECAEDRO (DUALE DELL'ICOSAEDRO)
	else if (p == 4 || p == 5)  
	{
		Polyhedron mesh = CreazioneSolidoPlatonico(p); // definisco il solido platonico di partenza

		if (b != c) // TRIANGOLAZIONE DI CLASSE I
		{
			cout<<"Eseguiamo la triangolazione di classe I"<<endl;
			unsigned int t_value = b + c; // valore che mi indica in quante parti dividere ogni lato del triangolo
			Polyhedron mesh_triangolata1 = TriangolazioneClasse1(mesh, t_value);	
			Polyhedron mesh_dualizzata = Dualizzazione(mesh_triangolata1);
			Polyhedron mesh_proiettata = ProiezioneSullaSfera(mesh_dualizzata);
			poliedro_finale = mesh_proiettata;
		}
		else // TRIANGOLAZIONE DI CLASSE II
		{
			cout<<"Eseguiamo la triangolazione di classe II"<<endl;
			Polyhedron mesh_triangolata2 = TriangolazioneClasse2(mesh, b);
			Polyhedron mesh_dualizzata = Dualizzazione(mesh_triangolata2);
			Polyhedron mesh_proiettata = ProiezioneSullaSfera(mesh_dualizzata);
			poliedro_finale = mesh_proiettata;
		}
	}

	
	// Cerco il cammino minimo	
	if(argc == 7)
	{	
		// Controllo validità v1 e v2
		if (v1 >= poliedro_finale.NumCell0Ds || v2 >= poliedro_finale.NumCell0Ds) 
		{
			cerr<<"I valori di v1 e v2 non sono validi: devono appartenere all'intervallo [0, " << poliedro_finale.NumCell0Ds -1<< "]" <<endl;
			return 1;
		}

		vector<Gedim::UCDProperty<double>> points_properties;
		vector<Gedim::UCDProperty<double>> segments_properties;
		points_properties.reserve(1); // Proprietà: ShortPath
		segments_properties.reserve(1);

		// Imposto inizialmente i valori di tutti i vertici e lati a 0; poi imposteremo ad 1 quelli che appartengono al cammino minimo
		vector<double> proprieta_vertici(poliedro_finale.NumCell0Ds, 0.0);
		vector<double> proprieta_lati(poliedro_finale.NumCell1Ds, 0.0);

		// Riempio point_properties
		Gedim::UCDProperty<double> prop_punto;
		prop_punto.Label = "ShortPath";
		prop_punto.NumComponents = 1;
		prop_punto.Data = proprieta_vertici.data();
		points_properties.push_back(prop_punto);

		// Riempio segments_properties
		Gedim::UCDProperty<double> prop_lato;
		prop_lato.Label = "ShortPath";
		prop_lato.NumComponents = 1;
		prop_lato.Data = proprieta_lati.data();
		segments_properties.push_back(prop_lato);


		auto [cammino, distanza] = CalcoloCamminoMinimo(poliedro_finale, v1, v2);
		cout<<"Il cammino minimo che collega "<< v1 <<" e "<< v2 <<" ha "<<(cammino.size()-1)<< " lati "<< endl;
        cout << "Il cammino minimo che collega " << v1 << " e " << v2 << " ha lunghezza: " << distanza << endl;
	
		// Riempio prop_vert
		for(unsigned int id_v : cammino) proprieta_vertici[id_v] = 1.0; //imposta a 1.0 la proprietà dei vertici lungo il percorso

		// Riempio prop_edges
		for(unsigned int i=0; i < cammino.size() - 1; i++)
		{
			unsigned int corrente = cammino[i];
			unsigned int successivo = cammino[i+1];
			unordered_set<unsigned int> punti_vicini = {corrente, successivo};

			// Itero sugli spigoli e, se lo spigolo è formato da (corrente, successivo), imposto prop_edges[j] = 1
			for(unsigned int j=0; j < poliedro_finale.Cell1DsExtrema.cols(); j++)
			{
				Vector2i segmento = poliedro_finale.Cell1DsExtrema.col(j);
				unsigned int a = segmento[0]; 
				unsigned int b = segmento[1];
				if(punti_vicini.count(a) + punti_vicini.count(b) == 2)
				{
					proprieta_lati[j] = 1.0;
				}
			}
		}
		EsportazionePoliedro(poliedro_finale, points_properties, segments_properties);
	}

	else
	{
		EsportazionePoliedro(poliedro_finale);
	}

	cout << "I dati del solido geodetico sono stati salvati nei file: 'Cell0Ds.txt','Cell1Ds.txt','Cell2Ds.txt','Cell3Ds.txt' " << endl;
	
	return 0;

}


