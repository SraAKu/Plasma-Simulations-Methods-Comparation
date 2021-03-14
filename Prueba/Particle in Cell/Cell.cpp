//Cell class Implementation//

#include "Cell.h"

#include <iostream>
using namespace std;

#include <math.h>

Cell::Cell()
{
	//cout << "o";
}

Cell::~Cell()
{
	//cout << "Celda destruida" << endl;
}

void Cell::set( Node *c1, Node *c2, Node *c3, Node *c4,
	            double ex, double ey, 
	            double ggx, double ggy,
	            int cellI)
{
	Ex = ex;
	Ey = ey;

	gridx = ggx;
	gridy = ggy;

	cellIndex = cellI;

	//cout << "Celda: " << cellIndex << " con esquinas: " << endl;

	corners[0] = c1;
	corners[1] = c2;
	corners[2] = c3;
	corners[3] = c4;

	/*for (int u = 0; u < 8; u++) {
		cout << " " << corners[u]->nodeIndex;
	}/**/

	//cout << endl;

}

void Cell::weightCharge( double q, double x, double y)
{
	//cout << "Acceso a weightCharge exitoso" << endl;
	//cout << x << " " << y << " " << z << endl;	
	double VT = gridx * gridy;
	//cout << "VT calculado" << endl;

	for (int b = 0; b < 4; b++) { 
		//cout << "Ponderacion en: " << b << endl;
		//cout << (x - corners[b]->gx()) << " ";
		//cout << (y - corners[b]->gy()) << " ";
		//cout << (z - corners[b]->gz()) << endl;
		//cout << "Inicia V:" << endl;
		double V = fabs( (x - corners[b]->gx()) * (y - corners[b]->gy()) );
		//cout << "V calculado" << endl;
		corners[b]->setCharge( (1 - V/VT) * q );
	}

}

void Cell::setField()
{
	//cout << "Acceso a setField exitoso" << endl;

	//comprobacion potenciales//
	/*for (int c = 0; c < 8; c++) {
		cout << "Potencial corner " << c+1 << ": ";
		cout << corners[c]->gPhi() << endl; 
	}
	/**/


	//comprobacion potenciales//



	Ex = ( -1 * (corners[1]->gPhi() - corners[0]->gPhi()) / gridx ) + 
	     ( -1 * (corners[3]->gPhi() - corners[2]->gPhi()) / gridx );


	//cout << "Ex = " << Ex;

	Ey = ( -1 * (corners[2]->gPhi() - corners[0]->gPhi()) / gridy ) + 
		 ( -1 * (corners[3]->gPhi() - corners[1]->gPhi()) / gridy );

	//cout << " Ey = " << Ey;


}

double Cell::gEx() { return Ex; }
double Cell::gEy() { return Ey; }

double Cell::ggridx() { return gridx; }
double Cell::ggridy() { return gridy; }
