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
	cout << argc << endl;
	if (argc == 5) {
		p = stoi(argv[1]);
		q = atoi(argv[2]);
		b = atoi(argv[3]);
		c = atoi(argv[4]);
        cout << "ok" << p << q << b << c << endl;
    }
	else if (argc == 7) {
		p = atoi(argv[1]);
		q = atoi(argv[2]);
		b = atoi(argv[3]);
		c = atoi(argv[4]);
		v1 = atoi(argv[5]);
		v2 = atoi(argv[6]);
        cout << "ok" << p << q << b << c << v1 << v2 << endl;
    }
	else
 	{
		cerr << "Inserimento non valido. Inserire 4 o 6 numeri interi positivi.\n" <<endl;
		return 1;
	}

	// controllo le condizioni su p e q
	cout << p << " " << q << " "<< double (1.0/p) <<endl;
	if (p < 3 || q < 3 || double (1.0/p)+double (1.0/q) <= 0.5){
		cerr << "Errore: i parametri p e q non rispettano le condizioni dei poligoni platonici.\n";
    	return 1;
	}
	

	// controllo la disugualianza dei poligoni platonici


	// controllo delle classi di triangolazione valide
	
	if (b != c && b!= 0 && c!= 0){
		cerr << "Errore: classe non trattata.\n";
		return 1;
	}

	std::cout << "Input ricevuto: p=" << p << ", q=" << q << ", b=" << b << ", c=" << c << "\n";




	Polyhedron mesh;
	return 0;
}