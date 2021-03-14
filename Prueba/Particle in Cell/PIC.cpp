//Particle in cell Test//

#include "Net.h"
#include "Cell.h"
#include "Node.h"
#include "Particle.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>

//#include "SFML/Graphics.hpp"

#include <iostream>
#include <fstream>
//using namespace std;

#include <chrono>

struct registro {
	float step;
	int pp;
	double ppx,ppy,ppvx,ppvy;
};

typedef registro Registro;

/*void PollEvent(sf::RenderWindow* pTarget, bool* pIsPaused, sf::View* pSimView);                                                         //Call all polled events for the sf::window
void Render(sf::RenderWindow* pTarget, std::vector<Particle*> pBodies, sf::Color pObjColor);                                                //Render given body objects to given screen
void SetView(sf::View* pView, sf::RenderWindow* pTarget, float pViewWidth, float pViewHeight);/**/                                        //set the window to the simulation view
//void DrawNode(Node* pNode, sf::RenderWindow* pTarget);                                                                                  //Draw a node to the screen, and all of its children (recursive)

void update(Net *, std::vector<Particle*> particles, double);
void PopulateBodyVectorUniform(std::vector<Particle*> *pBodies, unsigned int pParticlesCount, 
							   float pWidth, float pHeight, 
							   float pMinMass, float pMaxMass);

float const SimWidth = 500;                  //Width and height of simulation, needs to be large, particles outside of this range will not be included in the octree
float const SimHeight = 500;

float const ViewWidth = 1920;                   //Width and height of view of the simulation for the screen. 
float const ViewHeight = 1080;

/*float zoom = 5;                                 //The current amount of zoom in or out the user has inputed in total
float MouseX = 0;
float MouseY = 0;
bool IsPaused = false;                          //Contains the state of weather the simulation is paused or not
sf::Color ObjColor(255, 255, 255, 128);         //the defult colour of the objects
sf::View SimulationView;                        
sf::RenderWindow window(sf::VideoMode(1920, 1080), "N-Body simulation");
/**/
void PopulateBodyVectorUniform(std::vector<Particle*> *pBodies, unsigned int pParticlesCount, 
                               float pWidth, float pHeight, 
                               float pMinMass, float pMaxMass)
{
    srand(0);

    //printf("Acceso a populate exitoso\n");

    for (unsigned int i = 0; i < pParticlesCount; i++)
    {
    	//printf("particula %d\n", i );
    	Particle* p_i = new Particle();
        float positionx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * pWidth + (SimWidth / 2 - pWidth / 2);                         
        float positiony = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * pHeight + (SimHeight / 2 - pHeight / 2);

        float mass = pMinMass; //+ static_cast <float> (rand() % static_cast <int> (pMaxMass - pMinMass));                     //random mass (int) in range (MinObjectMass, MaxObjectMass)

        p_i->set(positionx,positiony,0,0,mass,mass,0);

        pBodies->push_back(p_i);
        //printf("Particula generada\n");
    }

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

	std::vector<Particle*> particles;

	/*Particle* p1 = new Particle();
	p1->set(1,1,1);
	Particle* p2 = new Particle();
	p2->set(0.3,1,1);
	Particle* p3 = new Particle();
	p3->set(0.1,0.1,0.1);
	Particle* p4 = new Particle();
	p4->set(1.5,0.2,0.5);
	Particle* p5 = new Particle();
	p5->set(1.1,1.3,1.5);
	Particle* p6 = new Particle();
	p6->set(1.1,1.7,0.4);
	Particle* p7 = new Particle();
	p7->set(0.8,0.2,1.7);
	Particle* p8 = new Particle();
	p8->set(1.8,1.8,1.3);

	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	particles.push_back(p6);
	particles.push_back(p7);
	particles.push_back(p8);/**/



	Net red(SimHeight,SimWidth,50,50);

	double dt = 1;

	FILE *ptrF;
	ptrF = fopen("Simulacion PIC.dat","w");
	Registro line;
	line.step = 0;

	PopulateBodyVectorUniform(&particles, nump, SimWidth, SimHeight, 1, 2);

	//printf("Numero de particulas: %d\n",particles.size() );

	for (int p = 0; p < particles.size(); p++) {
		line.pp = p;
		line.ppx = particles[p]->gx();
		line.ppy = particles[p]->gy();
		line.ppvx = particles[p]->gvx();
		line.ppvy = particles[p]->gvy();

		fprintf(ptrF, "%.2f %d %.4f %.4f %.4f %.4f \n",
		 	    line.step, line.pp, line.ppx, line.ppy, line.ppvx, line.ppvy);
		/*printf("%.2f %d %.4f %.4f %.4f %.4f %.4f %.4f\n",
		 	    line.step, line.pp, line.ppx, line.ppy, line.ppz, line.ppvx, line.ppvy, line.ppvz);
		printf("%s\n","casita");/**/
	}

	//SetView(&SimulationView, &window, ViewWidth, ViewHeight);

	for (int t = 1; t < steps; t++ ) {

		/*window.isOpen();
		PollEvent(&window, &IsPaused, &SimulationView); //These will always be done
		Render(&window, particles, ObjColor);
		/**/
		red.interaction(particles, particles.size());

		//cout << "Interaccion finalizada" << endl;

		red.C_PhiAndField();

		update(&red, particles, dt);

		//printf("%d\n",t);

		line.step = t * dt;

		for (int p = 0; p < particles.size(); p++) {
			line.pp = p;
			line.ppx = particles[p]->gx();
			line.ppy = particles[p]->gy();
			line.ppvx = particles[p]->gvx();
			line.ppvy = particles[p]->gvy();

			fprintf(ptrF, "%.2f %d %.4f %.4f %.4f %.4f\n",
			 line.step, line.pp, line.ppx, line.ppy, line.ppvx, line.ppvy);
			/*printf("%.2f %d %.4f %.4f %.4f %.4f %.4f %.4f\n",
			 line.step, line.pp, line.ppx, line.ppy, line.ppz, line.ppvx, line.ppvy, line.ppvz );/**/
		}

	}

	for (int pI = 0; pI < particles.size(); pI++) {
		delete particles[pI];
	}/**/

	//delete &red;

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() <<std::endl;

	std::ofstream outfile;
	outfile.open("output PIC.txt", std::ios_base::app);
	outfile << nump << " " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << std::endl;

	return 0;

}

void update(Net *net, std::vector<Particle*> particles, double dt)
{

	for (int pp = 0; pp < particles.size(); pp++) {

		particles[pp]->setvx( particles[pp]->gvx() + 
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


//****************Funciones de graficacion****************//

/*void PollEvent(sf::RenderWindow* pTarget, bool* pIsPaused, sf::View* pSimView)
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
/**/

/*void Render(sf::RenderWindow* pTarget, std::vector<Particle*> pBodies, sf::Color pObjColor)
{
    pTarget->clear();

    sf::RectangleShape Temp;
    //Temp.setFillColor(pObjColor);

    for (unsigned int i = 0; i < pBodies.size(); i++)
    {       
        if (zoom > 1)
            Temp.setSize(sf::Vector2f(pBodies.at(i)->gm() * zoom, pBodies.at(i)->gm() * zoom));
        else
            Temp.setSize(sf::Vector2f(pBodies.at(i)->gm(), pBodies.at(i)->gm()));   

        /*float ForceCoiffecent = sqrt(pBodies.at(i)->accX * pBodies.at(i)->accX + pBodies.at(i)->accY * pBodies.at(i)->accY) * (40000 * _GRAV_CONST);

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
        Temp.setPosition(pBodies.at(i)->gx(), pBodies.at(i)->gy());
        pTarget->draw(Temp);
    }

    //DrawNode(&GlobalNode, pTarget);

    pTarget->display();
}
/**/


/*void DrawNode(Node* pNode, sf::RenderWindow* pTarget)
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
/**/

/*void SetView(sf::View* pView, sf::RenderWindow* pTarget, float pViewWidth, float pViewHeight)
{
    pView->reset(sf::FloatRect(SimWidth / 2 - pViewWidth / 2, SimHeight / 2 - pViewHeight / 2, pViewWidth, pViewHeight));
    pView->setViewport(sf::FloatRect(0.f, 0.f, 1.f, 1.f));
    pTarget->setView(*pView);
}
/**/