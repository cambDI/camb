#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include "indigo.h"
#include "indigo-renderer.h"
#include "indigo-inchi.h"
#include "smindigo.h"

using namespace std;

extern "C" {
    
//uint32_t reverse(uint32_t x)
//{
//    x = ((x >> 1) & 0x55555555u) | ((x & 0x55555555u) << 1);
//    x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
//    x = ((x >> 4) & 0x0f0f0f0fu) | ((x & 0x0f0f0f0fu) << 4);
//    x = ((x >> 8) & 0x00ff00ffu) | ((x & 0x00ff00ffu) << 8);
//    x = ((x >> 16) & 0xffffu) | ((x & 0xffffu) << 16);
//    return x;
//}
    
int countHeavyNeighbours(int atom) {
    int numHeavyNeighbours = 0;
    int neighbour;
    int neighbourIterator = indigoIterateNeighbors(atom);

    while(neighbour = indigoNext(neighbourIterator)) {
        if(indigoAtomicNumber(neighbour) != 1) {
            numHeavyNeighbours++;
        }
        indigoFree(neighbour);  
    }
    indigoFree(neighbourIterator);
    
    return numHeavyNeighbours;
}

int countHydrogenNeighbours(int atom) {
    int numHydrogenNeighbours = 0;
    int neighbour;
    int neighbourIterator = indigoIterateNeighbors(atom);

    while(neighbour = indigoNext(neighbourIterator)) {
        if(indigoAtomicNumber(neighbour) == 1) {
            numHydrogenNeighbours++;
        }
        indigoFree(neighbour);  
    }
    indigoFree(neighbourIterator);
    
    return numHydrogenNeighbours;
}
    
void R_Fingerprints(char **structures_file, int *structureNumberIn, char **filename, int *useNameAsTitleInt) {
    int structure, fileIter;
    int structureNumber = structureNumberIn[0];
    printf("structureNumber: %d\n", structureNumber);
    bool useNameAsTitle = (*useNameAsTitleInt!=0);
    fileIter = indigoIterateSDFile(*structures_file);
    int i = 0;   
    while ((structure = indigoNext(fileIter))) {
        i++;
        if(i==structureNumber) {
            //indigoFoldHydrogens(structure);
            printf("filename: %s\n\n", *filename);
            if(useNameAsTitle) {
                renderMolecule(structure, *filename);
            }
            else {
                renderMolecule(structure, "", *filename, -1, "");
            }
            
            int atom, atomIter;
            atomIter = indigoIterateAtoms(structure);
            while(atom = indigoNext(atomIter)) {
                if(indigoAtomicNumber(atom) == 1) {
                    indigoFree(atom);
                    continue;
                }
                
                int atomicNumber = indigoAtomicNumber(atom);
                Rprintf("Atomic Number: %d\n\n", atomicNumber);
                
                int heavyNeighbours = countHeavyNeighbours(atom);
                Rprintf("Heavy Neighbours: %d\n\n", heavyNeighbours);
                
                int hydrogenNeighbours = countHydrogenNeighbours(atom);
                Rprintf("Hydrogen Neighbours: %d\n\n", hydrogenNeighbours);
                
                int charge = 0;
                indigoGetCharge(atom, &charge);
                Rprintf("Charge: %d\n\n", charge);
                
                //int isotope = indigoIsotope(atom);
                //Rprintf("Isotope: %d\n\n", isotope);
                
                
                indigoFree(atom);  
            }
            indigoFree(atomIter);
            
            
            break;
        }
        indigoFree(structure);
    }
    indigoFree(fileIter);
}

} // extern "C"
