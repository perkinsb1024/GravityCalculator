//
//  main.cpp
//  GravityCalculator
//
//  Created by Ben Perkins on 5/4/14.
//  Copyright (c) 2014 perkinsb1024. All rights reserved.
//

/* =====================
 TO DO:
 - Implement scaling (distance per pixel, velocity units, etc.)
 - Better data format
    - Easier to read (parse)
    - Compressed
    - FITS? http://en.wikipedia.org/wiki/FITS
 - Implement collisions
 - Implement materials
    - Elasticity
    - Density (?)
    - Drag coefficient (?)
 - Implement atmosphere (?)
 - Drop Python, implement OpenGL or similar
 
 DONE: 
 - Size independent of mass - 5/11/14
 ==================== */

#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>

const int STEPS = 1500; // Frames
const int RESOLUTION = 100; // Calculations per frame
const double G = 1;//6.67384 * pow(10, -11);
const double DELTA_T = 0.5 / RESOLUTION;
int CURRENT_PARTICLE = 0;

inline double random(double min, double max){ double f = (double)rand() / RAND_MAX; return min + f * (max - min); }

class _3vec {
private:
public:
    double x;
    double y;
    double z;
    
    _3vec();
    _3vec(double, double, double);
    
    // Assignment and equality
    inline void operator=(_3vec a) {x=a.x; y=a.y; z=a.z;};
    inline bool operator==(_3vec a) {return (x == a.x && y == a.y && z == a.z);};
    // Vector manipulations
    inline _3vec operator+(_3vec a) {return _3vec(x + a.x, y + a.y, z + a.z);};
    inline _3vec operator-(_3vec a) {return _3vec(x - a.x, y - a.y, z - a.z);};
    inline _3vec operator*(_3vec a) {return _3vec(x * a.x, y * a.y, z * a.z);};
    inline _3vec operator/(_3vec a) {return _3vec(x / a.x, y / a.y, z / a.z);};
    inline void operator+=(_3vec a) {*this = *this + a;};
    inline void operator-=(_3vec a) {*this = *this - a;};
    inline void operator*=(_3vec a) {*this = *this * a;};
    inline void operator/=(_3vec a) {*this = *this / a;};
    // Scalar manipulations
    inline _3vec operator+(double a) {return _3vec(x + a, y + a, z + a);};
    inline _3vec operator-(double a) {return _3vec(x - a, y - a, z - a);};
    inline _3vec operator*(double a) {return _3vec(x * a, y * a, z * a);};
    inline _3vec operator/(double a) {return _3vec(x / a, y / a, z / a);};
    inline void operator+=(double a) {*this = *this + a;};
    inline void operator-=(double a) {*this = *this - a;};
    inline void operator*=(double a) {*this = *this * a;};
    inline void operator/=(double a) {*this = *this / a;};
};

_3vec::_3vec(){};
_3vec::_3vec(double _x, double _y, double _z) { x = _x; y = _y; z = _z; }

class Force : public _3vec {
public:
    using _3vec::_3vec;
    Force(){};
    Force(_3vec a) {x=a.x; y=a.y; z=a.z;};
    inline void operator=(_3vec a) {x=a.x; y=a.y; z=a.z;};
};
class Acceleration : public _3vec {
public:
    using _3vec::_3vec;
    Acceleration(){};
    Acceleration(_3vec a) {x=a.x; y=a.y; z=a.z;};
    inline void operator=(_3vec a) {x=a.x; y=a.y; z=a.z;};
};
class Velocity : public _3vec {
public:
    using _3vec::_3vec;
    Velocity(){};
    Velocity(_3vec a) {x=a.x; y=a.y; z=a.z;};
    inline void operator=(_3vec a) {x=a.x; y=a.y; z=a.z;};
};
class Position : public _3vec {
public:
    using _3vec::_3vec;
    Position(){};
    Position(_3vec a) {x=a.x; y=a.y; z=a.z;};
    inline void operator=(_3vec a) {x=a.x; y=a.y; z=a.z;};
    double distance(Position p);
};

double Position::distance(Position p) {
    return sqrt(pow((x-p.x), 2) +
                pow((y-p.y), 2) +
                pow((z-p.z), 2));
}


class Particle {
private:
    double mass;
    double size;
    Velocity velocity;
    Position position;
    int id;
public:
    Particle(double _mass, Velocity _velocity, Position _position);
    Particle(double _mass, double _size, Velocity _velocity, Position _position, int _r, int _g, int _b);
    inline bool operator==(Particle a);
    void setMass(double _mass);
    double getMass();
    void setSize(double _size);
    double getSize();
    void setVelocity(Velocity _velocity);
    Velocity getVelocity();
    void setPosition(Position _position);
    Position getPosition();
    int getId();
    bool isTouching(Particle p);
    void collide(Particle p);
    
    // todo: Fix
    int r;
    int g;
    int b;
    
    void print();
    void print(std::string prefix);
};

Particle::Particle(double _mass, Velocity _velocity, Position _position) {
    id = CURRENT_PARTICLE++;
    mass = _mass;
    size = _mass;
    velocity = _velocity;
    position = _position;
    r = 255;
    g = 255;
    b = 255;
}

Particle::Particle(double _mass, double _size, Velocity _velocity, Position _position, int _r, int _g, int _b) {
    id = CURRENT_PARTICLE++;
    mass = _mass;
    size = _size;
    velocity = _velocity;
    position = _position;
    r = _r;
    g = _g;
    b = _b;
}

inline bool Particle::operator==(Particle a) {
    return id == a.getId();
}

double Particle::getMass() {
    return mass;
}

void Particle::setMass(double _mass) {
    mass = _mass;
}

void Particle::setSize(double _size) {
    size = _size;
}

double Particle::getSize() {
    return size;
}

void Particle::setVelocity(Velocity _velocity) {
    velocity = _velocity;
}

Velocity Particle::getVelocity() {
    return velocity;
}

void Particle::setPosition(Position _position) {
    position = _position;
}

Position Particle::getPosition() {
    return position;
}

int Particle::getId() {
    return id;
}

bool Particle::isTouching(Particle p) {
    return position.distance(p.getPosition()) <= size + p.getSize();
}

void Particle::collide(Particle p) {
    // Note: A collision between particles 'a' and 'b' must be called for particles:
        // a.collide(b);
        // b.collide(a);
    // a.collide(b) will not modify particle 'b' in any way
    
    // Note to self, this will not work unless particles store their pre-collision velocities separately from their post-collision velocities
}

void Particle::print() {
    print("");
}

void Particle::print(std::string prefix) {
    printf("%sMass: %f\n%sPosition: {%f, %f, %f}\n%sVelocity:{%f, %f, %f}\n", prefix.c_str(), mass, prefix.c_str(), position.x, position.y, position.z, prefix.c_str(), velocity.x, velocity.y, velocity.z);
}

// todo: Declare rest of functions:

void translateToCenterOfMass(std::vector<Particle> &universe);
void centerOnPosition(std::vector<Particle> &universe, Position p);
Particle randomParticle();
Particle randomParticle(double, double, double, double, Velocity, Velocity, Position, Position, int, int, int);

Force F_g(Particle a, Particle b) {
    if(a == b) return Force(0, 0, 0);
    double x, y, z;
    double total, distance, theta, phi;
    Position aPos = a.getPosition();
    Position bPos = b.getPosition();
    x = aPos.x - bPos.x;
    y = aPos.y - bPos.y;
    z = aPos.z - bPos.z;
    distance = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    total = G*a.getMass()*b.getMass()/sqrt(distance);
    theta = asin(abs(y)/distance);
    phi = asin(abs(z)/distance);
    
    return Force ((x < 0 ? 1 : -1) * total * cos(theta),
                    (y < 0 ? 1 : -1) * total * sin(theta),
                    (z < 0 ? 1 : -1) * total * sin(phi));
}

Force computeForce(std::vector<Particle> universe, Particle compare) {
    Force force (0, 0, 0);
    for(int i = 0; i < universe.size(); i++) {
        force+=F_g(compare, universe.at(i));
    }
    return force;
}

Velocity incrementVelocity(Velocity v, Acceleration a) {
    return Velocity( v.x += a.x * DELTA_T, v.y += a.y * DELTA_T, v.z += a.z * DELTA_T);
}

Position incrementPosition(Position p, Velocity v) {
    return Position( p.x += v.x * DELTA_T, p.y += v.y * DELTA_T, p.z += v.z * DELTA_T);
}

Acceleration accelerationFromForce(Force f, double m) {
    return Acceleration( f.x / m, f.y / m, f.z / m );
}

void updateVelocities(std::vector<Particle> &universe) {
    for(int i = 0; i < universe.size(); i++) {
        Force f = Force(0, 0, 0);
        Particle p = universe[i];
        double m = p.getMass();
        for(int j = 0; j < universe.size(); j++){
            f += F_g(p, universe.at(j));
        }
        p.setVelocity(incrementVelocity(p.getVelocity(), accelerationFromForce(f, m)));
        universe[i] = p;
    }
}

void updatePositions(std::vector<Particle> &universe) {
    for(int i = 0; i < universe.size(); i++) {
        Particle p = universe[i];
        p.setPosition(incrementPosition(p.getPosition(), p.getVelocity()));
        universe[i] = p;
    }
}

Position getCenterOfMass(std::vector<Particle> universe) {
    Position com = Position(0, 0, 0);
    double mass = 0;
    for(int i = 0; i < universe.size(); i++) {
        Particle p = universe.at(i);
        mass += p.getMass();
        com += p.getPosition() * p.getMass();
    }
    return com / mass;
}

bool orderByDepth(Particle a, Particle b) {
    return a.getPosition().z < b.getPosition().z;
}

void printUniverse(std::vector<Particle> universe, int time) {
    std::cout << "=======Time " << time << "=======" << "\n";
    std::sort(universe.begin(), universe.end(), orderByDepth);
    for(int i = 0; i < universe.size(); i++) {
        Particle p = universe.at(i);
        //std::cout << "Particle " << i << "\n";
        //p.print(" ");
        Position pos = p.getPosition();
        //printf("%f,%f,%f,%f,%d,%d,%d\n", p.getSize(), pos.x, pos.y, pos.z, p.r, p.g, p.b);
        printf("%d,%d,%d,%d,%d,%d,%d\n", (int)p.getSize(), (int)pos.x, (int)pos.y, (int)pos.z, p.r, p.g, p.b);
    }
    
    //Position com = getCenterOfMass(universe);
    //printf("CoM: %f,%f,%f\n", com.x, com.y, com.z);
    //printf("CoM: %d,%d,%d\n", (int)com.x, (int)com.y, (int)com.z);
}

void translateToCenterOfMass(std::vector<Particle> &universe) {
    centerOnPosition(universe, getCenterOfMass(universe));
}

void centerOnPosition(std::vector<Particle> &universe, Position p) {
    for(int j = 0; j < universe.size(); j++) {
        universe.at(j).setPosition(universe.at(j).getPosition() - p);
    }
}

Particle randomParitcle() {
    return randomParticle(0.01, 10, 1, 10, Velocity(-1, -1, -1), Velocity(1, 1, 1), Position(-500, -300, -100), Position(500, 300, 100), -1, -1, -1);
}

Particle randomParticle(double minMass, double maxMass,
                        double minSize, double maxSize,
                        Velocity minVelocity, Velocity maxVelocity,
                        Position minPosition, Position maxPosition,
                        int r, int g, int b) {
    
    if(r < 0 || r > 255) { r = (int)random(0, 255); }
    if(g < 0 || g > 255) { g = (int)random(0, 255); }
    if(b < 0 || b > 255) { b = (int)random(0, 255); }
    
    return Particle(
                     random(minMass, maxMass), // Mass
                     random(minSize, maxSize), // Size
                     Velocity(
                                  random(minVelocity.x, maxVelocity.x), // v_x
                                  random(minVelocity.y, maxVelocity.y), // v_y
                                  random(minVelocity.z, maxVelocity.z)  // v_z
                              ),
                     Position(
                                  random(minPosition.x, maxPosition.x), // s_x
                                  random(minPosition.y, maxPosition.y), // s_y
                                  random(minPosition.z, maxPosition.z)  // s_z
                              ),
                     r, // red
                     g, // green
                     b  // blue
                    );
}

int main(int argc, const char * argv[])
{
    srand (uint(time(NULL))); // Seed randomness
    
    std::vector<Particle> universe;
    /*
    universe.push_back(Particle(40, Velocity(0, 0, 0), Position(-200, 200, 100), 255, 0, 0));
    universe.push_back(Particle(30, Velocity(0, 0, 0), Position(0, 0, 20), 0, 255, 0));
    universe.push_back(Particle(10, Velocity(-0.1, 0, 0), Position(5, 1, 1), 255, 255, 0));
    universe.push_back(Particle(10, Velocity(0, 0.1, 0), Position(1, 2, 3), 0, 0, 255));
    */
    
    // Sun, Earth, Moon... ish
    universe.push_back(Particle(330, 100, Velocity(0, 0, 0), Position(0, 0, 0), 255, 200, 0));
    universe.push_back(Particle(0.001, 15, Velocity(0, -75, 0), Position(300, 0, 0), 0, 50, 255));
    universe.push_back(Particle(0.000012, 5, Velocity(0, -70, 2), Position(320, 0, 0), 255, 255, 255));
    
    Particle *earth = &universe[1];
    
    /*
    for(int i = 0; i < 50; i++) {
        universe.push_back(randomParticle);
    }
     */
    
    for(int i = 0; i < STEPS; i++) {
        for(int j = 0; j < RESOLUTION; j++) {
            updateVelocities(universe);
            updatePositions(universe);
            centerOnPosition(universe, earth->getPosition());
        }
        printUniverse(universe, i);
    }
    
    return 0;
}

