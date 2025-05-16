#include <iostream>

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <cmath>
#include <algorithm>

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
	if (p == 3)
	{
		auto [verts, faces] = getSolidData(q);
		Polyhedron mesh = buildPlatonicSolid(p, q, b, c);
	
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

		if (b != c) //TRINAGOLAZIONE DI CLASSE I
		{
			cout<<"triangolazione di classe I"<<endl;
			unsigned int t_value = b + c; //valore che mi indica in quante parti dividere ogni lato del triangolo
		}

		else //TRINAGOLAZIONE DI CLASSE II
		{
			cout<<"triangolazione di classe II"<<endl;
		}

	}

	else if (p == 4)  //CASO CUBO (DUALE DELL'OTTAEDRO)
		{cout<<"p = 4"<<endl;}

	else if (p == 5)  //CASO DODECAEDRO (DUALE DELL'ICOSAEDRO)
		{cout<<"p = 5"<<endl;}


	return 0;

}
		
	

