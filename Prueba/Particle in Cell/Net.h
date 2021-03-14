//Net class Definition//

#ifndef NET_H
#define NET_H

#include <vector>

#include "Node.h"
#include "Cell.h"

#include "Particle.h"

class Net {
public:
	Net( double = 1, double = 1, 
		 double = 1, double = 1);
	~Net();
	void interaction( std::vector<Particle*> particles, int N);
	void C_PhiAndField();
	int gpx();
	int gpy();
	double gCEx( int );
	double gCEy( int );
	int fV( int );

private:
	double px, py;
	int Nx, Ny;
	int nx, ny;

	std::vector<Node*> nodes;
	std::vector<Cell*> cells;
	/*Node nodes[ 9261 ];
	Cell cells[ 8000 ];/**/

};

#endif
