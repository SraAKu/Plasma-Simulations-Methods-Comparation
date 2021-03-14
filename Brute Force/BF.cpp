//Brute Force//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>

#include <chrono>

#define _PI 3.14159265      //Pi, used for calculations and rounded to 8 decimal places. 
#define _GRAV_CONST 0.1 //6.674e-11     //the gravitational constant. This is the timestep between each frame. Lower for slower but more accurate simulations

struct Body
{
    float posX, posY;           //position x and y
    float velX, velY;           //velocity x and y
    double accX, accY;      //force acting on object since last frame
    float mass;                 //mass of object
};

struct registro {
    float step;
    int pp;
    double ppx,ppy,ppvx,ppvy;
};

typedef registro Registro;

void CalculateForce(Body* bi, Body* bj, float pSoftener);
void ResuelveNodo(std::vector<Body*> Bodies, float pSoftener);  
Body* CreateBody(float px, float py, float pmass, float pvx = 0, float pvy = 0);                                                        //return a pointer to new body object defined on the heap with given paramiters                                                                            //Calculate force exerted on eachother between two bodies
void PopulateBodyVectorUniform(std::vector<Body*> *pBodies, unsigned int pParticlesCount, float pWidth, float pHeight, float pMinMass, float pMaxMass);
void UpdateBodies(std::vector<Body*> pBodies);  
void ResetForces(std::vector<Body*> pBodies);                                                                                        //Calculate velocity chance from the bodies force exerted since last update, update position based on velocity, reset force to 0
void DeleteBodies(std::vector<Body*> pBodies);                                                                                          //Deletes objects pointed to by given vector

std::vector<Body*> Bodies; 

float const SimWidth = 500;                  //Width and height of simulation, needs to be large, particles outside of this range will not be included in the octree
float const SimHeight = 500;
float const Softener = 10; 
bool IsPaused = false;                      

float dt = 1;

inline void CalculateForce(Body* bi, Body* bj, float pSoftener)  //bi is being attracted to bj. 15 flops of calculation
{
    //std::vector from i to j
    float vectorx = bj->posX - bi->posX;
    float vectory = bj->posY - bi->posY;

    //c^2 = a^2 + b^2 + softener^2
    float distSqr = vectorx * vectorx + vectory * vectory + pSoftener * pSoftener;

    // ivnDistCube = 1/distSqr^(3/2)
    float distSixth = distSqr * distSqr * distSqr;
    double invDistCube = 1.0f / (sqrt(distSixth));


    double acc = (bj->mass * invDistCube * _GRAV_CONST);

    bi->accX += vectorx * acc;
    bi->accY += vectory * acc;
}

inline void ResuelveNodo(std::vector<Body*> Bodies, float pSoftener)
{
    for (unsigned int i = 0; i < Bodies.size(); i++) {
        for (unsigned int j = i; j < Bodies.size(); j++) {
            CalculateForce(Bodies.at(i),Bodies.at(j), Softener);
        }
    }
}

Body* CreateBody(float px, float py, float pmass, float pvx, float pvy)
{
    Body* Temp = new Body;

    Temp->posX = px;
    Temp->posY = py;
    Temp->mass = pmass;
    Temp->velX = pvx;
    Temp->velY = pvy;
    Temp->accX = Temp->accY = 0;

    return Temp;
}

//Distribucion de las particulas de manera uniforme//
void PopulateBodyVectorUniform(std::vector<Body*> *pBodies, unsigned int pParticlesCount, 
                               float pWidth, float pHeight, 
                               float pMinMass, float pMaxMass)
{
    srand(0);

    for (unsigned int i = 0; i < pParticlesCount; i++)
    {
        float positionx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * pWidth + (SimWidth / 2 - pWidth / 2);                         
        float positiony = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * pHeight + (SimHeight / 2 - pHeight / 2);

        float mass = pMinMass; //+ static_cast <float> (rand() % static_cast <int> (pMaxMass - pMinMass));                     //random mass (int) in range (MinObjectMass, MaxObjectMass)

        pBodies->push_back(CreateBody(positionx, positiony, mass));
    }
}

void UpdateBodies(std::vector<Body*> pBodies)
{
    for (unsigned int i = 0; i < pBodies.size(); i++)
    {
        while (pBodies[i]->posX > SimWidth)
        {  
            pBodies[i]->posX -= SimWidth;
        }

        while (pBodies[i]->posX < 0)
        {
            pBodies[i]->posX += SimWidth;
        }

        while (pBodies[i]->posY > SimHeight)
        {
            pBodies[i]->posY -= SimHeight;
        }

        while (pBodies[i]->posY < 0)
        {
            pBodies[i]->posY += SimHeight;
        }

        pBodies.at(i)->velX += pBodies.at(i)->accX;     //f = ma => force = mass * accelleration. Therefor
        pBodies.at(i)->velY += pBodies.at(i)->accY;     //a = f/m => accelleration = force / mass

        pBodies.at(i)->posX += pBodies.at(i)->velX;
        pBodies.at(i)->posY += pBodies.at(i)->velY;
    }
}

void ResetForces(std::vector<Body*> pBodies)
{
    for (unsigned int i = 0; i < pBodies.size(); i++)
    {
        pBodies[i]->accX = 0;
        pBodies[i]->accY = 0;
    }
}

void DeleteBodies(std::vector<Body*> pBodies)
{
    for (unsigned int i = 0; i < pBodies.size(); i++)
    {
        delete pBodies[i];
    }

    pBodies.clear();
}


int main()
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    int i;
    int steps;
    float nump;
    std::ifstream inFile;
    
    inFile.open("input.txt");
    
    for (int i = 0; i < 2; i++) {
        inFile >> steps;
        inFile >> nump;
    }

    /*PopulateBodyVectorDisk(&Bodies, NumParticles, 
                           SimWidth, SimHeight, 
                           DiskRadiusMax, DiskRadiusMin, 
                           ObjectMassMin, ObjectMassMax, 
                           GalaticCenterMass);/**/

    PopulateBodyVectorUniform(&Bodies, nump, 
                               SimWidth, SimHeight, 
                               1, 2);/**/

    /*Bodies.push_back(CreateBody(0, 0, 10, 100, 100));
    Bodies.push_back(CreateBody(0, 0, 10, 100, -100));
    Bodies.push_back(CreateBody(0, 0, 10, -100, 100));*/

    FILE *ptrF;
    ptrF = fopen("Simulacion BF.dat","w");
    Registro line;
    line.step = 0;

    for (int p = 0; p < Bodies.size(); p++) {
        line.pp = p;
        line.ppx = Bodies[p]->posX;
        line.ppy = Bodies[p]->posY;
        line.ppvx = Bodies[p]->velX;
        line.ppvy = Bodies[p]->velY;

        fprintf(ptrF, "%.2f %d %.4f %.4f %.4f %.4f \n",
                line.step, line.pp, line.ppx, line.ppy, line.ppvx, line.ppvy);
        /*printf("%.2f %d %.4f %.4f %.4f %.4f %.4f %.4f\n",
                line.step, line.pp, line.ppx, line.ppy, line.ppz, line.ppvx, line.ppvy, line.ppvz);
        printf("%s\n","casita");/**/
    }



    //SetView(&SimulationView, &window, ViewWidth, ViewHeight);

    for (int i = 0; i < steps; i++)
    {
        //PollEvent(&window, &IsPaused, &SimulationView); //These will always be done

        if (!IsPaused)  //These will not if the simulation is paused
        {       
            //AttractToCenter(Bodies, SimWidth, SimHeight, GalaticCenterMass);
            UpdateBodies(Bodies);
            ResetForces(Bodies);
            ResuelveNodo(Bodies, Softener);
            //GlobalNode.Reset();
            //GlobalNode.SetParam(Bodies, SimWidth, SimHeight);
            //OctreeBodyAttraction();
        }   

        //Render(&window, Bodies, ObjColor);

        line.step = i * dt;

        for (int p = 0; p < Bodies.size(); p++) {
            line.pp = p;
            line.ppx = Bodies[p]->posX;
            line.ppy = Bodies[p]->posY;
            line.ppvx = Bodies[p]->velX;
            line.ppvy = Bodies[p]->velY;

            fprintf(ptrF, "%.2f %d %.4f %.4f %.4f %.4f\n",
             line.step, line.pp, line.ppx, line.ppy, line.ppvx, line.ppvy);
            /*printf("%.2f %d %.4f %.4f %.4f %.4f %.4f %.4f\n",
             line.step, line.pp, line.ppx, line.ppy, line.ppz, line.ppvx, line.ppvy, line.ppvz );/**/
        }
    }

    DeleteBodies(Bodies);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() <<std::endl;

    std::ofstream outfile;
    outfile.open("output BF.txt", std::ios_base::app);
    outfile << nump << " " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << std::endl;


    return 0;
}