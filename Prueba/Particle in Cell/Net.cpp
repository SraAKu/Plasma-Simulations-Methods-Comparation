//Net class Implementation//

#include "Net.h"

#include <iostream>
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fasttransforms.h"
using namespace std;

#define _PI 3.14159265 

Net::Net( double x, double y, 
	      double cx, double cy)
{
	px = x;
	py = y;
	Nx = cx * 1;
	Ny = cy * 1;
	nx = Nx + 1;
	ny = Ny + 1;
	double ggx = px / cx;
	double ggy = py / cy;

	int i = 0;
	int j = 0;
	int k = 0;

	//cout << "\nInicia construccion de nodos" << endl;

	for ( int d = 0; d < (nx * ny); d++) {

		double xx = i * ggx;
		double yy = j * ggy;
		Node* n = new Node();
		n->set( xx, yy,
			    0 , 0, 1,
				nx * (k * ny + j) + i);

		nodes.push_back(n);

		/*nodes[d]->set( xx, yy, zz,
					  0 , 0, 1,
					  nx * (k * ny + j) + i);/**/

		i += 1;

		if ( i >= nx ) {
			i = 0;
			j += 1;
		};

	};

	//cout << "inicia construccion ultimos dos nodos" << endl;

	int size = nodes.size();
	Node* m = new Node();
	m->set(-1.0, -1.0, 
			0, 0, 1,
			size);

	Node* l = new Node();
	l->set(-1.0, -1.0,
			0, 0, 1,
			size + 1);

	nodes.push_back(m);
	nodes.push_back(l);

	/*nodes[nx*ny*nz].set( -1.0, -1.0, -1.0,
						 0, 0, 1,
						 9259);/**/

	//cout << "Inicia construccion celdas" << endl;

	int control = 0;

	for ( int c = 0; c < (Nx * Ny); c++) {

		if ( (c >= Nx) && (c % Nx == 0)) {
			control++;
		}

		if ( (c >= (Nx*Ny)) && ((c % (Nx*Ny)) == 0 )) {
			control = control + nx; 
		}

		Cell* celli = new Cell();

		celli->set( nodes[fV(control)], nodes[fV(control+1)], nodes[fV(control+nx)], nodes[fV(control+nx+1)],
		             0.0, 0.0, 
		             ggx, ggy,
		             c);

		cells.push_back(celli);

		control++;
	};

	//cout << "---------- Red construida ----------" << endl;

}

Net::~Net()
{
	//cout << "Destruyendo Red" << endl;

	for (int n = 0; n <= nx * ny + 1; n++ ) {
		delete nodes[n];
	};

	for (int N = 0; N < Nx * Ny; N++ ) {
		delete cells[N];
	};/**/

	//cout << "Red destruida" << endl;
}

void Net::interaction(std::vector<Particle*> particles, int N)
{
	//cout << "Inicia interaccion" << endl;

	for (int p = 0; p < N; p++) {

		int i = Nx * particles[p]->gx() / px;
		int j = Ny * particles[p]->gy() / py;
		int k = 0;
		int index = Nx * ( k * Ny + j ) + i;
		//cout << "indice particula " << p << ": " << i << " " << j << " " << k << "= " << index << endl;

		particles[p]->seti(index);

		double Q = particles[p]->gq();
		double X = particles[p]->gx();
		double Y = particles[p]->gy();

		cells[index]->weightCharge(Q, X, Y);

	};

	//cout << "Ponderacion de carga finalizada" << endl;

}

void Net::C_PhiAndField()
{
	double Phis[nx*ny];
	double KG_pi = 0.1 / _PI; //2.12440018e-11;
	int dim = nx*ny;

	for (int nn = 0; nn < nx*ny; nn++) {
		Phis[nn] = nodes[nn]->gQ();
	};

	//printf("Matriz construida\n");

	alglib::real_1d_array tempA;
	tempA.setlength(nx);
	alglib::complex_1d_array I;
	I.setlength(dim);
	alglib::complex_1d_array tempI;
	tempI.setlength(nx);
	alglib::complex_1d_array II;
	II.setlength(dim);
	alglib::complex_1d_array tempII;
	tempII.setlength(nx);
	alglib::complex_1d_array F;
	F.setlength(dim);
	int cont;
	int n;

	//printf("Arreglos definidos\n");

	//FFT en Filas//
	cont = 0;
	int fila = 0;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {

			n = nx * (0 * ny + j) + i;

			tempA[cont] = Phis[n];
			cont++;

			if ( i == nx-1 ) {
				cont = 0;
				fftr1d(tempA, tempI);

				for (int p = 0; p < nx; p++) {
					I[fila*nx+p] = tempI[p];
				};

				fila++;

			};	

		};

	};
	

	//printf("Primera transformada exitosa\n" );

	//FFT en columnas//
	cont = 0;
	int columna = 0;

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {

			n = nx * (0 * ny + j) + i;

			tempI[cont] = I[n];
			cont++;

			if ( j == ny-1 ) {
				cont = 0;
				fftc1d(tempI);

				for (int p = 0; p < ny; p++) {
					II[columna*ny+p] = tempI[p];
				};

				columna++;

			};

		};
	};

	//printf("Segunda transformada exitosa\n" );

	cont = 0;
	int final = 0;

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {

			n = nx * (0 * ny + j) + i;

			tempI[cont] = I[n];
			cont++;

			if ( j == ny-1 ) {
				cont = 0;

				for (int p = 0; p < ny; p++) {
					F[final*ny+p] = tempI[p];
				};

				final++;

			};

		};
	};

	//printf("*<--->*Re-ordenacion*<--->*\n" );

	double V2 = px*py*px*py;
	double lx2 = px*px;
	double ly2 = py*py;
	double alfa;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {

			n = nx * (0 * ny + j) + i;

			if ( n == 0 ) {
				alfa = 1;
			}
			else {
				alfa = - KG_pi * V2 / ((ly2 * i*i ) + (lx2 * j*j));
			}

			F[n] = alfa * F[n];

			};
		};

	//printf("Transformacion Phi exitosa\n" );

	//***********************INVERSA*******************//

	//IFFT en Filas//
	cont = 0;
	fila = 0;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {

			n = nx * (0 * ny + j) + i;

			tempI[cont] = F[n];
			cont++;

			if ( i == nx-1 ) {
				cont = 0;
				fftc1dinv(tempI);

				for (int p = 0; p < nx; p++) {
					I[fila*nx+p] = tempI[p];
				};

				fila++;

			};

		};
	};

	//printf("Primera transformada inversa exitosa\n" );

	//IFFT en columnas//
	cont = 0;
	columna = 0;

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {

			n = nx * (0 * ny + j) + i;

			tempII[cont] = I[n];
			cont++;

			if ( j == ny-1 ) {
				cont = 0;
				fftc1dinv(tempII);

				for (int p = 0; p < ny; p++) {
					II[columna*ny+p] = tempII[p];

				};

				columna++;

			};

		};
	};

	//printf("Segunda transformada inversa exitosa\n" );

	cont = 0;
	final = 0;

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {

			n = nx * (0 * ny + j) + i;

			tempI[cont] = I[n];
			cont++;

			if ( j == ny-1 ) {
				cont = 0;

				for (int p = 0; p < ny; p++) {
					nodes[final*ny+p]->setPhi(tempI[p].x);
				};

				final++;

			};

		};
	};

	for ( int ii = 0; ii < nx * ny; ii++) {
		nodes[ii]->reCharge();
	};

	//cout << "Cargas re-establecidas" << endl;

	for (int II = 0; II < Nx * Ny; II++) {
		//cout << "campo celda " << II << ": " << endl;
		cells[II]->setField();
	};

	//cout << "Campo establecido" << endl;

}

int Net::fV( int i )
{
	if ( ( i < 0 ) || (i >= nx * ny) ) {
		return nx * ny;
	};

	return i;

}

int Net::gpx() { return px; }
int Net::gpy() { return py; }

double Net::gCEx( int u ) { return cells[u]->gEx(); }
double Net::gCEy( int v ) { return cells[v]->gEy(); }
