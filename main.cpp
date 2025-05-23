#include <iostream>
#include "Polyhedron.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

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
	if (argc == 5) 
	{
		p = stoi(argv[1]);  //commenta stoi
		q = stoi(argv[2]);
		b = stoi(argv[3]);
		c = stoi(argv[4]);
		cout << "Input ricevuto: p=" << p << ", q=" << q << ", b=" << b << ", c=" << c << "\n";
    }
	else if (argc == 7) 
	{
		p = stoi(argv[1]);
		q = stoi(argv[2]);
		b = stoi(argv[3]);
		c = stoi(argv[4]);
		v1 = stoi(argv[5]);
		v2 = stoi(argv[6]);
		cout << "Input ricevuto: p=" << p << ", q=" << q << ", b=" << b << ", c=" << c << ", v1=" << v1 << ", v2=" << v2 << "\n";
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
	if (p == 3)
	{	
		Polyhedron mesh = buildPlatonicSolid(q); //definiamo il solido platonico di partenza
		// export del poliedro di partenza
		Gedim::UCDUtilities utilities;
			{
				utilities.ExportPoints("./Cell0Ds.inp",
									mesh.Cell0DsCoordinates);
			}

			{
				utilities.ExportSegments("./Cell1Ds.inp",
										mesh.Cell0DsCoordinates,
										mesh.Cell1DsExtrema);
			}

		if (b != c) //TRIANGOLAZIONE DI CLASSE I
		{
			cout<<"triangolazione di classe I"<<endl;
			unsigned int t_value = b + c; //valore che mi indica in quante parti dividere ogni lato del triangolo
			Polyhedron mesh_triangolata1 = triangulateClass1(mesh, t_value); //triangolazione di classe 1 con parametro t_value

			// export del poliedro con triangolzione di classe 1
			Gedim::UCDUtilities utilities;
				{
					utilities.ExportPoints("./Cell0Ds_T1.inp",
										mesh_triangolata1.Cell0DsCoordinates);
				}

				{
					utilities.ExportSegments("./Cell1Ds_T1.inp",
											mesh_triangolata1.Cell0DsCoordinates,
											mesh_triangolata1.Cell1DsExtrema);
				}
		}

		else //TRIANGOLAZIONE DI CLASSE II
		{
			cout<<"triangolazione di classe II"<<endl;
			/*
			Polyhedron mesh_triangolata2 = triangulateClass2(mesh, b);

			// export del poliedro con triangolzione di classe 1
			Gedim::UCDUtilities utilities;
				{
					utilities.ExportPoints("./Cell0Ds_T2.inp",
										mesh_triangolata2.Cell0DsCoordinates);
				}

				{
					utilities.ExportSegments("./Cell1Ds_T2.inp",
											mesh_triangolata2.Cell0DsCoordinates,
											mesh_triangolata2.Cell1DsExtrema);
				}
			*/
		}
	}

	//CASO CUBO (DUALE DELL'OTTAEDRO) e CASO DODECAEDRO (DUALE DELL'ICOSAEDRO)
	else if (p == 4 || p == 5)  
		{
			cout<<"p = "<< p <<endl;
			Polyhedron mesh = buildPlatonicSolid(p); //definiamo il solido platonico di partenza

			if (b != c) //TRIANGOLAZIONE DI CLASSE I
			{
				cout<<"triangolazione di classe I"<<endl;
				unsigned int t_value = b + c; //valore che mi indica in quante parti dividere ogni lato del triangolo
				Polyhedron mesh_triangolata1 = triangulateClass1(mesh, t_value); //triangolazione di classe 1 con parametro t_value
					
				Polyhedron mesh_dualizzata = Dualize(mesh_triangolata1);

				// export del poliedro con dualizzazione
				Gedim::UCDUtilities utilities;
					{
						utilities.ExportPoints("./Cell0Ds_D1.inp",
											mesh_dualizzata.Cell0DsCoordinates);
					}

					{
						utilities.ExportSegments("./Cell1Ds_D1.inp",
												mesh_dualizzata.Cell0DsCoordinates,
												mesh_dualizzata.Cell1DsExtrema);
					}
			}

			else //TRIANGOLAZIONE DI CLASSE II
			{
				cout<<"triangolazione di classe II"<<endl;
			}
		
		}

	return 0;
}


		
	

