//Cell class Definition//

#ifndef CELL_H
#define CELL_H

#include "Node.h"

class Cell {
public:
	Cell();
	~Cell();
	void set( Node *, Node *, Node *, Node *,
	          double = 0.0, double = 0.0,
	          double = 1.0, double = 1.0,
	          int = 0);
	void weightCharge(double = 0.0, double = 0.0, 
			          double = 0.0);
	void setField();
	double gEx();
	double gEy();
	double ggridx();
	double ggridy();

	int cellIndex;

private:
	Node* corners[4];
	double Ex, Ey;
	double gridx, gridy;

};

#endif
