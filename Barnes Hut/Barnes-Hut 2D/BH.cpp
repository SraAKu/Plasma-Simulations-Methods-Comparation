#include "SFML/Graphics.hpp"
#include "Node.h"
#include <ctime>
#include <math.h>
#include <vector>

#include <iostream>
#include <fstream>

#define _PI 3.14159265      //Pi, used for calculations and rounded to 8 decimal places. 
#define _GRAV_CONST 0.1     //the gravitational constant. This is the timestep between each frame. Lower for slower but more accurate simulations

struct registro {
    float step;
    int pp;
    double ppx,ppy,ppvx,ppvy;
};

typedef registro Registro;

void BodyAttraction(std::vector<Body*> pBodies, float pSoftener);                                                                       //Attracts each body to each other body in the given vector of pointers to body objects
void CalculateForceNode(Body* bi, Node* bj, float pSoftener);                                                                           //Calculate force exerted on body from node
void CalculateForce(Body* bi, Body* bj, float pSoftener);
void ResuelveNodo(Node* N, float pSoftener);                                                                               //Calculate force exerted on eachother between two bodies
Body* CreateBody(float px, float py, float pmass, float pvx = 0, float pvy = 0);                                                        //return a pointer to new body object defined on the heap with given paramiters
void DeleteBodies(std::vector<Body*> pBodies);                                                                                          //Deletes objects pointed to by given vector
void PollEvent(sf::RenderWindow* pTarget, bool* pIsPaused, sf::View* pSimView);                                                         //Call all polled events for the sf::window
void PopulateBodyVectorDisk(std::vector<Body*> *pBodies, unsigned int pParticlesCount, float pWidth, float pHeight, float pMaxDist, float pMinDist, float pMinMass, float pMaxMass, float pGalaticCenterMass = 0);  //populate given vector with bodies with given paramiters in a disk formation
void PopulateBodyVectorUniform(std::vector<Body*> *pBodies, unsigned int pParticlesCount, float pWidth, float pHeight, float pMinMass, float pMaxMass);
void Render(sf::RenderWindow* pTarget, std::vector<Body*> pBodies, sf::Color pObjColor);                                                //Render given body objects to given screen
void SetView(sf::View* pView, sf::RenderWindow* pTarget, float pViewWidth, float pViewHeight);                                          //set the window to the simulation view
void UpdateBodies(std::vector<Body*> pBodies);                                                                                          //Calculate velocity chance from the bodies force exerted since last update, update position based on velocity, reset force to 0
void DrawNode(Node* pNode, sf::RenderWindow* pTarget);                                                                                  //Draw a node to the screen, and all of its children (recursive)
void CheckNode(Node* pNode, Body* pBody);                                                                                               //Checks if a node is sufficently far away for force calculation, if not recureses on nodes children
void OctreeBodyAttraction();                                                                                                            //Using a calculated oct-tree, calculate the force exerted on each object
void AttractToCenter(std::vector<Body*> pBodies, float width, float height, float centerMass);                                                                                                                  //Attract each particle to the center of the simulation
void ResetForces(std::vector<Body*> pBodies);
void RepellFromCenter(std::vector<Body*> pBodies, float width, float height, float centerMass);

float const DiskRadiusMax = 300;              //Max and min distances objects will be from the galatic center
float const DiskRadiusMin = 50;
float const GalaticCenterMass = 1000000;        //The mass of the very large object simulating a black hole at the center of a galixy;
float const ObjectMassMax = 2;                  //The max and min mass of the objects in the galixy
float const ObjectMassMin = 1;
float const SimWidth = 500;                  //Width and height of simulation, needs to be large, particles outside of this range will not be included in the octree
float const SimHeight = 500;
float const ViewWidth = 1920;                   //Width and height of view of the simulation for the screen. 
float const ViewHeight = 1080;
float const Softener = 10;                      //A softener used for the force calculations, 10 is a good amount 
float dt = 1;     //Number of particles in simtulation, currently 2^15                                
double const _NODE_THRESHOLD = 100;             //Threshold for node calculations   

float zoom = 5;                                 //The current amount of zoom in or out the user has inputed in total
float MouseX = 0;
float MouseY = 0;

std::vector<Body*> Bodies;                      //Container of all Bodies in simulations
Node GlobalNode;

bool IsPaused = false;                          //Contains the state of weather the simulation is paused or not
sf::Color ObjColor(255, 255, 255, 128);         //the defult colour of the objects
sf::View SimulationView;                        
sf::RenderWindow window(sf::VideoMode(1920, 1080),"N-Body simulation");

//Distribucion de las particulas en disco alrededor de una masa central//
void PopulateBodyVectorDisk(std::vector<Body*> *pBodies, unsigned int pParticlesCount, 
                            float pWidth, float pHeight, 
                            float pMaxDist, float pMinDist, 
                            float pMinMass, float pMaxMass, 
                            float pGalaticCenterMass)
{
    srand(0);

    for (unsigned int i = 0; i < pParticlesCount; i++)
    {
        float angle = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2 * static_cast <float> (_PI))));    //sets angle to random float range (0, 2 pi) 

        float distanceCoefficent = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float distance = pMinDist + ((pMaxDist - pMinDist) * (distanceCoefficent * distanceCoefficent));                    //Distance point will be from the galatic center, between MinDiskRadius and MaxDiskRadius

        float positionx = cos(angle) * distance + (pWidth / 2);                                                             //set positionx and positiony to be the point you get when you go in the direction of 'angle' till you have traveled 'distance' 
        float positiony = sin(angle) * distance + (pHeight / 2);

        float orbitalVelocity = sqrt((pGalaticCenterMass * static_cast <float> (_GRAV_CONST)) / distance);                  //Calculate the orbital velocity required to orbit the galatic centre   

        float velocityx = (sin(angle) * orbitalVelocity);
        float velocityy = (-cos(angle) * orbitalVelocity);

        float mass = pMinMass + static_cast <float> (rand() % static_cast <int> (pMaxMass - pMinMass));                     //random mass (int) in range (MinObjectMass, MaxObjectMass)

        pBodies->push_back(CreateBody(positionx, positiony, mass, velocityx, velocityy));                                   
    }
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

        float mass = pMinMass + static_cast <float> (rand() % static_cast <int> (pMaxMass - pMinMass));                     //random mass (int) in range (MinObjectMass, MaxObjectMass)

        pBodies->push_back(CreateBody(positionx, positiony, mass));
    }
}

//Genera cuerpo y retorna un aputador al mismo//
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

//Calcula fuerza generada en cuerpo bi por cuerpo bj//
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

//Resuelve nodo//
inline void ResuelveNodo(Node* N, float pSoftener)
{
    for (unsigned int i = 0; i < (N->Bodies).size(); i++) {
        for (unsigned int j = i; j < (N->Bodies).size(); j++) {
            CalculateForce((N->Bodies).at(i),(N->Bodies).at(j), Softener);
        }
    }
}

//Calcula fuerza generada en cuerpo bi por centro de masa de nodo bj//
inline void CalculateForceNode(Body* bi, Node* bj, float pSoftener)  //bi is being attracted to bj. 15 flops of calculation
{
    //vector from the body to the center of mass
    float vectorx = bj->CenterOfMassx - bi->posX;
    float vectory = bj->CenterOfMassy - bi->posY;

    //c^2 = a^2 + b^2 + softener^2
    float distSqr = vectorx * vectorx + vectory * vectory + pSoftener * pSoftener;

    // ivnDistCube = 1/distSqr^(3/2)
    float distSixth = distSqr * distSqr * distSqr;
    double invDistCube = 1.0f / (sqrt(distSixth));


    double acc = (bj->TotalMass * invDistCube * _GRAV_CONST);

    bi->accX += vectorx * acc;
    bi->accY += vectory * acc;
}

//Calcula fuerza generada por el cuerpo central sobre todas las particulas//
void AttractToCenter(std::vector<Body*> pBodies, float width, float height, float centerMass)
{
    Body* Temp = CreateBody(width / 2, height / 2, centerMass); //Create a body at the center of the simulation

    for (unsigned int i = 0; i < pBodies.size(); i++)
    {
        CalculateForce(pBodies[i], Temp, Softener);
    }

    delete Temp;
}

//Actualiza velocidades y posiciones de particulas//
void UpdateBodies(std::vector<Body*> pBodies)
{
    for (unsigned int i = 0; i < pBodies.size(); i++)
    {
        if ((pBodies[i]->posX > SimWidth && pBodies[i]->velX > 0) || (pBodies[i]->posX < 0 && pBodies[i]->velX < 0))
        {
            pBodies[i]->velX = -pBodies[i]->velX;
        }

        if ((pBodies[i]->posY > SimHeight && pBodies[i]->velY > 0) || (pBodies[i]->posY < 0 && pBodies[i]->velY < 0))
        {
            pBodies[i]->velY = -pBodies[i]->velY;
        }

        pBodies.at(i)->velX += pBodies.at(i)->accX;     //f = ma => force = mass * accelleration. Therefor
        pBodies.at(i)->velY += pBodies.at(i)->accY;     //a = f/m => accelleration = force / mass

        pBodies.at(i)->posX += pBodies.at(i)->velX;
        pBodies.at(i)->posY += pBodies.at(i)->velY;
    }
}

/*
void UpdateBodies(std::vector<Body*> particles)
{

    for (int pp = 0; pp < particles.size(); pp++) {

        particles[pp]->velX += particles( particles[pp]->gvx() + 
                         net->gCEx(particles[pp]->gi()) * dt );

        particles[pp]->setvy( particles[pp]->gvy() + 
                         net->gCEy(particles[pp]->gi()) * dt );

        double x = particles[pp]->gx() + particles[pp]->gvx() * dt;
        particles[pp]->setx(x);

        while ( x >= net->gpx() ) {
            x = x - net->gpx();
            particles[pp]->setx(x);
        }

        while ( x < 0 ) {
            x = net->gpx() + x;
            particles[pp]->setx(x);
        }

        double y = particles[pp]->gy() + particles[pp]->gvy() * dt;
        particles[pp]->sety(y);

        while ( y >= net->gpy() ) {
            y = y - net->gpy();
            particles[pp]->sety(y);
        }

        while ( y < 0 ) {
            y = net->gpy() + y;
            particles[pp]->sety(y);
        }

    }

}
/**/

//Reinicia fuerzas//
void ResetForces(std::vector<Body*> pBodies)
{
    for (unsigned int i = 0; i < pBodies.size(); i++)
    {
        pBodies[i]->accX = 0;
        pBodies[i]->accY = 0;
    }
}

//Funcion recursiva que recorre el arbol calculando la fuerza sobre un cuerpo//
inline void CheckNode(Node* pNode, Body* pBody)
{
    if (pNode->Bodies.size() != 0)
    {
        float diffX = (pNode->CenterOfMassx - pBody->posX);
        float diffY = (pNode->CenterOfMassy - pBody->posY);

        float distance = sqrt((diffX) * (diffX) + (diffY) * (diffY));   //Distance from the node to the object                          

        if ((pNode->width / distance) < (_NODE_THRESHOLD) || (pNode->HasChildren == false))     //if sufficently far away or has no children (external node) group node and calculate force
        {   
            CalculateForceNode(pBody, pNode, Softener);

            if (!pNode->IsUsed) {
                ResuelveNodo(pNode,Softener);
                pNode->IsUsed = true;
            }
        }
        else                                                                                    //if not, repeat function with children
        {
            CheckNode(pNode->Child[0], pBody);
            CheckNode(pNode->Child[1], pBody);
            CheckNode(pNode->Child[2], pBody);
            CheckNode(pNode->Child[3], pBody);
        }
    }
} 

void OctreeBodyAttraction()
{
    for (unsigned int i = 0; i < Bodies.size(); i++)
    {   
        CheckNode(&GlobalNode, Bodies[i]);
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

//Funcion de calculo de fuerzas entre particulas de un mismo nodo//
void BodyAttraction(std::vector<Body*> pBodies, float pSoftener)
{
    for (unsigned int i = 0; i < pBodies.size(); i++)
    {
        for (unsigned int j = 0; j < pBodies.size(); j++)
        {
            CalculateForce(pBodies.at(i), pBodies.at(j), pSoftener); //for each body in pBodies: each other body in pBodies: Calculate attractive force exerted on the first body from the second one
        }
    }
}

void RepellFromCenter(std::vector<Body*> pBodies, float width, float height, float centerMass)
{
    Body* Temp = CreateBody(width / 2, height / 2, centerMass); //Create a body at the center of the simulation

    for (unsigned int i = 0; i < pBodies.size(); i++)
    {
        float vectorx = Temp->posX - pBodies[i]->posX;
        float vectory = Temp->posY - pBodies[i]->posY;

        float distSqr = vectorx * vectorx + vectory * vectory;

        double Dist = (sqrt(distSqr));

        double acc = (pBodies[i]->mass * Dist * _GRAV_CONST * 0.0001);

        pBodies[i]->accX -= vectorx * acc;
        pBodies[i]->accY -= vectory * acc;
    }

    delete Temp;
}

int main()
{

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
                               ObjectMassMin, ObjectMassMax);/**/

    /*Bodies.push_back(CreateBody(0, 0, 10, 100, 100));
    Bodies.push_back(CreateBody(0, 0, 10, 100, -100));
    Bodies.push_back(CreateBody(0, 0, 10, -100, 100));*/

    FILE *ptrF;
    ptrF = fopen("Simulacion.dat","w");
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



    SetView(&SimulationView, &window, ViewWidth, ViewHeight);

    for (int i = 0; i < steps; i++)
    {
        //PollEvent(&window, &IsPaused, &SimulationView); //These will always be done

        if (!IsPaused)  //These will not if the simulation is paused
        {       
            //AttractToCenter(Bodies, SimWidth, SimHeight, GalaticCenterMass);
            UpdateBodies(Bodies);
            ResetForces(Bodies);
            GlobalNode.Reset();
            GlobalNode.SetParam(Bodies, SimWidth, SimHeight);
            OctreeBodyAttraction();
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
}

//****************Funciones de graficacion****************//

void PollEvent(sf::RenderWindow* pTarget, bool* pIsPaused, sf::View* pSimView)
{
    sf::Event event;

    while (pTarget->pollEvent(event))
    {
        if (event.type == sf::Event::Closed)
            pTarget->close();
        if (event.type == sf::Event::KeyPressed)
        {
            if (event.key.code == sf::Keyboard::Space)
                *pIsPaused = !*pIsPaused;                   //toggle what is pointed to by IsPaused
        }
        if (event.type == sf::Event::MouseWheelScrolled)
        {
            zoom *= 1 + (static_cast <float> (-event.mouseWheelScroll.delta) / 10); //for each notch down, -10%, for each notch up, +10%
            pSimView->zoom(1 + (static_cast <float> (-event.mouseWheelScroll.delta) / 10));
        }
    }

    if (sf::Mouse::getPosition().x > (1920 - 20))
        SimulationView.move(2 * zoom, 0);
    if (sf::Mouse::getPosition().x < (0 + 20))
        SimulationView.move(-2 * zoom, 0);
    if (sf::Mouse::getPosition().y > (1080 - 20))
        SimulationView.move(0, 2 * zoom);
    if (sf::Mouse::getPosition().y < (0 + 20))
        SimulationView.move(0, -2 * zoom);

    pTarget->setView(*pSimView);
}

void Render(sf::RenderWindow* pTarget, std::vector<Body*> pBodies, sf::Color pObjColor)
{
    pTarget->clear();

    sf::RectangleShape Temp;
    //Temp.setFillColor(pObjColor);

    for (unsigned int i = 0; i < pBodies.size(); i++)
    {       
        if (zoom > 1)
            Temp.setSize(sf::Vector2f(pBodies.at(i)->mass * zoom, pBodies.at(i)->mass * zoom));
        else
            Temp.setSize(sf::Vector2f(pBodies.at(i)->mass, pBodies.at(i)->mass));   

        float ForceCoiffecent = sqrt(pBodies.at(i)->accX * pBodies.at(i)->accX + pBodies.at(i)->accY * pBodies.at(i)->accY) * (40000 * _GRAV_CONST);

        if (ForceCoiffecent > 1)
            ForceCoiffecent = 1;

        float Red, Green, Blue;

        Blue = 1 - (ForceCoiffecent);

        if (ForceCoiffecent < 0.2)
            Red = (ForceCoiffecent) * 5;
        else
            Red = 1;

        if (ForceCoiffecent < 0.5)
            Green = (ForceCoiffecent) * 2;
        else
            Green = 1;

        Temp.setFillColor(sf::Color(Red * 255, Green * 255, Blue * 255, 128));
        Temp.setPosition(pBodies.at(i)->posX, pBodies.at(i)->posY);
        pTarget->draw(Temp);
    }

    //DrawNode(&GlobalNode, pTarget);

    pTarget->display();
}

void DrawNode(Node* pNode, sf::RenderWindow* pTarget)
{
        sf::RectangleShape Temp;
        Temp.setFillColor(sf::Color(0, 0, 0, 0));
        Temp.setOutlineThickness(zoom);
        Temp.setOutlineColor(sf::Color(0, 255, 0, 16));
        Temp.setPosition(pNode->posX, pNode->posY);
        Temp.setSize(sf::Vector2f(pNode->width, pNode->height));

        pTarget->draw(Temp);
    if (pNode->HasChildren) //recercivly draw all children
    {
        DrawNode(pNode->Child[0], pTarget);
        DrawNode(pNode->Child[1], pTarget);
        DrawNode(pNode->Child[2], pTarget);
        DrawNode(pNode->Child[3], pTarget);
    }
}

void SetView(sf::View* pView, sf::RenderWindow* pTarget, float pViewWidth, float pViewHeight)
{
    pView->reset(sf::FloatRect(SimWidth / 2 - pViewWidth / 2, SimHeight / 2 - pViewHeight / 2, pViewWidth, pViewHeight));
    pView->setViewport(sf::FloatRect(0.f, 0.f, 1.f, 1.f));
    pTarget->setView(*pView);
}
