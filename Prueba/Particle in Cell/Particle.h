//Particle class Definition//

#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {

public:
	Particle();
	~Particle();
	void set( double = 0.0, double = 0.0, 
		  double = 0.0, double = 0.0,
		  double = 0.01, double = 1e-4, int = 0);
	void setx( double );
	void sety( double );
	double gx();
	double gy();
	void setvx( double );
	void setvy( double );
	double gvx();
	double gvy();
	double gm();
	double gq();
	void seti( int );
	int gi();

private:
	double x, y;
	double vx, vy;
	double m, q;
	int i;

};

#endif
