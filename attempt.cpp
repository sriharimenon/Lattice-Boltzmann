#include "palabos2D.h"
#include "palabos2D.hh"
#include<iostream>
#include<iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

void writeGifs(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice),
                               imSize, imSize );
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
 
    const plint maxIter = 10000; // Iterate during 1000 steps.
    const plint nx = 100;       // Choice of lattice dimensions.
    const plint ny = 20;
    const T omega = 1.;        // Choice of the relaxation parameter
    const T density = 1.;
 
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
 
    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
           nx, ny, new BGKdynamics<T,DESCRIPTOR>(omega) );
 
    Array<T, 2> zeroVelocity(0,0);
    Array<T, 2> boundaryVelocity(0.01,0);
 
    //Inlet
    boundaryCondition->setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0,0, 1,ny-2) );
    //Outlet
    boundaryCondition->setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(nx-1,nx-1, 1,ny-2), boundary::outflow );
 
    // Velocity boundary condition on bottom wall.
    setBoundaryVelocity(lattice, Box2D(0, nx-1, 0, 0), Array<T,2>(0.,0.) );

    // Velocity boundary condition on top wall.
    setBoundaryVelocity(lattice, Box2D(0, nx-1, ny-1, ny-1), Array<T,2>(0.,0.) );
 
    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            boundaryVelocity );
 
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), density, zeroVelocity);
 
    
    lattice.initialize();

// Main loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT) {
        if (iT%100==0) {
            pcout << "Saving Gif at time step " << iT << endl;
            writeGifs(lattice, iT);
        }
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
