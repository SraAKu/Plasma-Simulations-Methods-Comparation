//Node class implementation//

#include "Node.h"

#include <iostream>
using namespace std;

Node::Node()
{
	//cout << ".";
}

Node::~Node()
{
	//cout << "Nodo destruido" << endl;
}

void Node::set( double xx, double yy,
	            double QQ, double PPhi, double eps,
	            int nodeI)
{
	Q = QQ;
	Phi = PPhi;
	x = xx;
	y = yy;
	epsilon = eps;

	nodeIndex = nodeI;
}

void Node::setCharge(double q) { Q += q; }
void Node::reCharge() { Q = 0.0; }
double Node::gQ() { return Q; }

void Node::setPhi(double p) { Phi += p; }
double Node::gPhi() { return Phi; }

double Node::gx() { return x; }
double Node::gy() { return y; }
