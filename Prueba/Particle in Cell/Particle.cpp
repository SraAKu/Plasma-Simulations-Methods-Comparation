//Particle class Implementation//

#include "Particle.h"

#include <iostream>
using namespace std;

Particle::Particle()
{
	//cout << "Particula construida" << endl;
}

Particle::~Particle()
{
	//cout << "Particula destruida" << endl;
}

void Particle::set( double xx, double yy,
					double vvx, double vvy, 
				    double mm, double qq, int a)
{
	x = xx;
	y = yy;
	vx = vvx;
	vy = vvy;
	m = mm;
	q = qq;
	i = a;
}

void Particle::setx( double tx ) { x = tx; }
void Particle::sety( double ty ) { y = ty; }

double Particle::gx() { return x; }
double Particle::gy() { return y; }

void Particle::setvx( double tvx ) { vx = tvx; }
void Particle::setvy( double tvy ) { vy = tvy; }

double Particle::gvx() { return vx; }
double Particle::gvy() { return vy; }

double Particle::gm() { return m; }
double Particle::gq() { return q; }

void Particle::seti( int a ) { i = a; }
int Particle::gi() { return i; }
