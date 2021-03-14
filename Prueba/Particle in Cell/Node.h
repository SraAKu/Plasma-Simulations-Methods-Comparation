//Node class Definition//

#ifndef NODE_H
#define NODE_H

class Node {
public:
	Node();
	~Node();
	void set( double = -2.0, double = -2.0,
		  double = 0.0, double = 0.0, double = 1.0,
		  int = -1);
	void setCharge( double q );
	void reCharge();
	double gQ();
	void setPhi( double p );
	double gPhi();
	double gx();
	double gy();

	int nodeIndex;

private:
	double Q, Phi;
	double x, y;
	double epsilon;

};

#endif
