#include <iostream>
#include "Polyhedron.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedronLibrary;

int main()
{
	unsigned int p;
	unsigned int q;
	unsigned int b;
	unsigned int c;
	unsigned int v1;
	unsigned int v2;
	cout << "Inserisci la quadrupla/sestupla del poliedro nel formato p q b c v1 v2 (lasciare vuoto nel caso non si vogliano inserire i valori di v1 e v2):\n";
	//cin >> p >> q >> b >> c;

	// controllo che il tipo dei valori inseriti sia int
	if (!(std::cin >> p >> q >> b >> c >> v1 >> v2)) {
        cerr << "Errore: inserire 4 numeri interi positivi.\n";
        return 1;
    }

	// controllo le condizioni su p e q
	if (p < 3 || q <3){
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