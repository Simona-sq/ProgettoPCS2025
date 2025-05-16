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
	if (argc == 5) {
		p = stoi(argv[1]);  //commenta stoi
		q = stoi(argv[2]);
		b = stoi(argv[3]);
		c = stoi(argv[4]);
		cout << "Input ricevuto: p=" << p << ", q=" << q << ", b=" << b << ", c=" << c << "\n";
    }
	else if (argc == 7) {
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
	if (p < 3 || q < 3 || double (1.0/p)+double (1.0/q) <= 0.5){
		cerr << "Errore: i parametri p e q non rispettano le condizioni dei poligoni platonici.\n";
    	return 1;
	}

	// controllo delle classi di triangolazione valide
	if (b != c && b!= 0 && c!= 0){
		cerr << "Errore: classe non trattata.\n";
		return 1;
	}


	auto [verts, faces] = getSolidData(q);
	Polyhedron mesh = buildPlatonicSolid(p, q, b, c);

	// stampa vertici
	std::cout << "Vertici:\n";
	for (size_t i = 0; i < verts.size(); ++i) 
	{
		const auto& v = verts[i];
		std::cout << "v" << i << ": (" 
				<< v(0) << ", " 
				<< v(1) << ", " 
				<< v(2) << ")\n";
	}

	// stampa facce
	std::cout << "Facce:\n";
	for (size_t i = 0; i < faces.size(); ++i) {
		std::cout << "f" << i << ": ";
		for (unsigned int vid : faces[i]) {
			std::cout << vid << " ";
		}
		std::cout << "\n";
	}
        
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


	//stampa id
	for (unsigned int i = 0; i < mesh.NumCell0Ds; ++i) 
	{
        std::cout << "ID " << i << " : "
                  << mesh.Cell0DsCoordinates(0, i) << " "
                  << mesh.Cell0DsCoordinates(1, i) << " "
                  << mesh.Cell0DsCoordinates(2, i) << "\n";
    }
		
	return 0;
}