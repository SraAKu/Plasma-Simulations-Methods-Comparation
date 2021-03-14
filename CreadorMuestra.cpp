#include <stdio.h>
#include <stdlib.h>
#include <fstream>

int const nump = 100;

float const SimWidth = 500;                  //Width and height of simulation, needs to be large, particles outside of this range will not be included in the octree
float const SimHeight = 500;

int main()
{
	FILE *ptrF;
	ptrF = fopen("Simulacion.dat","w");

	srand(0);

	for (unsigned int p = 0; p < nump; p++) {


		float positionx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * SimWidth;                         
    	float positiony = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * SimHeight;

    	fprintf(ptrF, "%.4f %.4f %.4f %.4f \n",
		 	    		positionx, positiony, 0.0, 0.0);

    }



}