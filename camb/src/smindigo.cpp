#include "./smindigo.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "indigo.h"
#include "indigo-renderer.h"
#include "indigo-inchi.h"

void renderMolecule(int structure, const char* folder, const char* file, int number, const char* title) {
    // render the molecule to file
    indigoLayout(structure);
    indigoSetOption("render-output-format", "png");
    indigoSetOption("render-comment", title);
    indigoSetOption("render-comment-position", "top");
    indigoSetOptionXY("render-image-size", 400, 500);
    indigoSetOptionColor("render-background-color", 1.0, 1.0, 1.0);
    char str[256];
    if(number != -1) {
        sprintf(str, "%s/%d%s%s", folder, number, "_", file);
    }
    else if(strcmp(folder,"")!=0)
    {
        sprintf(str, "%s/%s", folder, file);
    }
    else {
        sprintf(str, "%s", file);
    }
    indigoRenderToFile(structure, str);
}

void renderMolecule(int structure, const char* file) {
    const char* folder = "";
    int number = -1;
    const char* title = indigoName(structure);
    renderMolecule(structure, folder, file, -1, title);
}

void renderPair(int s1, int s2, const char* folder, const char* file, int number, const char* title) {
    int collection = indigoCreateArray();
    indigoArrayAdd(collection, s1);
    indigoArrayAdd(collection, s2);
    
    indigoSetOption("render-output-format", "png");
    indigoSetOption("render-comment", title);
    indigoSetOption("render-comment-position", "top");
    indigoSetOptionXY("render-image-size", 800, 1000);
    indigoSetOptionColor("render-background-color", 1.0, 1.0, 1.0);
    char str[256];
    if(number != -1) {
        sprintf(str, "./%s/%d%s%s%s", folder, number, "_", file, ".png");
    }
    else
    {
        sprintf(str, "%s%s", file, ".png");
    }
    indigoRenderGridToFile(collection, NULL, 2, str);
    indigoFree(collection);
}

// looks for CO- SO- and PO- instances, neutralises the O
// returns the number of hits
int neutraliseCSPnegO(int structure) {
    int query, matcher, matchIter, match, queryAtomIter, queryAtom; // Indigo objects that require freeing
    
    query = indigoLoadSmartsFromString("[O-$(*[#6,#14,#15])]");
    matcher = indigoSubstructureMatcher(structure, 0);
    matchIter = indigoIterateMatches(matcher, query);
    int numMatches = 0;
    if (matchIter != -1) {
        while (match = indigoNext(matchIter)) {
            queryAtomIter = indigoIterateAtoms(query);
            while (queryAtom = indigoNext(queryAtomIter))
            {
                int structureAtom = indigoMapAtom(match, queryAtom);
                indigoSetCharge(structureAtom, 0);
                indigoFree(queryAtom);
            }
            numMatches++;
            indigoFree(queryAtomIter);
            indigoFree(match);
        }
    }
    indigoFree(matchIter);
    indigoFree(matcher);
    indigoFree(query);
    return numMatches;
}

// looks for NH+ instances, neutralises the N
// returns the number of hits
int neutraliseNposH(int structure) {
    int query, matcher, matchIter, match, queryAtomIter, queryAtom; // Indigo objects that require freeing
    
    query = indigoLoadSmartsFromString("[#7+$(*[H])]"); // #7 means any nitrogen, aliphatic or aromatic 
    matcher = indigoSubstructureMatcher(structure, 0);
    matchIter = indigoIterateMatches(matcher, query);
    int numMatches = 0;
    if (matchIter != -1) {
        while (match = indigoNext(matchIter)) {
            queryAtomIter = indigoIterateAtoms(query);
            while (queryAtom = indigoNext(queryAtomIter))
            {
                int structureAtom = indigoMapAtom(match, queryAtom);
                indigoSetCharge(structureAtom, 0);
                indigoFree(queryAtom);
            }
            numMatches++;
            indigoFree(queryAtomIter);
            indigoFree(match);
        }
    }
    indigoFree(matchIter);
    indigoFree(matcher);
    indigoFree(query);
    return numMatches;
}

// returns the largest component within the structure (you need to free this returned object later)
// doesn't free structure
int pickLargestSubstructure(int structure) {
    int componentIter, component, largestComponent; // Indigo objects that require freeing

    componentIter = indigoIterateComponents(structure);
    int maxAtoms = 0;
    if (componentIter != -1) {
        while (component = indigoNext(componentIter)) {
            // if only one component
            if (indigoIndex(component) == 0 && !indigoHasNext(componentIter))
            {
                indigoFree(componentIter);
                return indigoClone(structure);
            }
            
            int currentComponent = indigoClone(component);
            int numAtoms = indigoCountAtoms(currentComponent);
            //Rprintf("Number of atoms: %d\n\n", numAtoms);
            if (numAtoms>maxAtoms) {
                maxAtoms = numAtoms;
                largestComponent = indigoClone(currentComponent);  
            }
            indigoFree(currentComponent);
            indigoFree(component);
        }
    }
    indigoFree(componentIter);
    return largestComponent;
}

bool isOrganic(int structure) {
    int atom, atomIter;
    bool structureOrganic = true;
    
    atomIter = indigoIterateAtoms(structure);
    while(atom = indigoNext(atomIter)) {
        
        bool atomOrganic = indigoAtomicNumber(atom) == 1 ||
                           indigoAtomicNumber(atom) == 6 ||
                           indigoAtomicNumber(atom) == 7 ||
                           indigoAtomicNumber(atom) == 8 ||
                           indigoAtomicNumber(atom) == 15 ||
                           indigoAtomicNumber(atom) == 16 ||
                           indigoAtomicNumber(atom) == 9 ||
                           indigoAtomicNumber(atom) == 17 ||
                           indigoAtomicNumber(atom) == 35 ||
                           indigoAtomicNumber(atom) == 53;
                           
        if(!atomOrganic) {
            structureOrganic = false;
            indigoFree(atom);
            break;
        }
        indigoFree(atom);  
    }
    indigoFree(atomIter);
    return structureOrganic;
}

bool isAromatic(int structure) {
    int bond, bondIter;
    bool structureAromatic = false;
    
    bondIter = indigoIterateBonds(structure);
    while(bond = indigoNext(bondIter)) {
        
        bool bondAromatic  = (indigoBondOrder(bond) == 4);
                           
        if(bondAromatic) {
            structureAromatic = true;
            indigoFree(bond);
            break;
        }
        indigoFree(bond);  
    }
    indigoFree(bondIter);
    return structureAromatic;
}

void resetIsotopes(int structure) {
    int atom, atomIter;
   
    atomIter = indigoIterateAtoms(structure);
    while(atom = indigoNext(atomIter)) {
        indigoResetIsotope(atom);
        indigoFree(atom);  
    }
    indigoFree(atomIter);
}

bool containsMoreThanX(int structure, int X, int atomicNumber) {
    int count = 0;
    int atom, atomIter;
   
    atomIter = indigoIterateAtoms(structure);
    while(atom = indigoNext(atomIter)) {
        if(indigoAtomicNumber(atom)==atomicNumber) {
            count++;
        }        
        indigoFree(atom);  
    }
    indigoFree(atomIter);
    return count > X;
}

