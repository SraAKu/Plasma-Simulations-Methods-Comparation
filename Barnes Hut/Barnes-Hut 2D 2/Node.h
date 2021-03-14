#pragma once
#include <vector>

struct Body
{
    float posX, posY;           //position x and y
    float velX, velY;           //velocity x and y
    double accX, accY;      //force acting on object since last frame
    float mass;                 //mass of object
};

class Node
{
public:
    Node();
    Node(unsigned int pdepth);
    ~Node();
    void GenerateChildren();
    void SetParam(std::vector<Body*> pBodies, float pwidth, float pheight, float px = 0, float py = 0);
    void Reset();

    std::vector<Body*> Bodies;
    std::vector<Node*> Child;
    bool HasChildren;

    float posX, posY;
    float width, height;
    float TotalMass;
    float CenterOfMassx;
    float CenterOfMassy;
    unsigned int Depth;
    bool IsUsed;            //For testing, delete this later
};
