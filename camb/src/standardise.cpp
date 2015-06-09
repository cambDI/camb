#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <boost/algorithm/string.hpp>
//#include <boost/lexical_cast.hpp>


#include "indigo.h"
#include "indigo-renderer.h"
#include "indigo-inchi.h"
#include "smindigo.h"

using namespace std;

string trimInchi(const char* inchi) {
    string trimmed(inchi);
    size_t found;
    
    found = trimmed.find("/t");
    if(found != string::npos) {
        trimmed = trimmed.substr(0, found);
    }
    
    found = trimmed.find("/b");
    if(found != string::npos) {
        trimmed = trimmed.substr(0, found);
    }
    return trimmed;
}

extern "C" {
    
    void R_drawMoleculeInSDF(char **structures_file, int *structureNumberIn, char **filename, int *useNameAsTitleInt) {
        int structure, fileIter;
        int structureNumber = structureNumberIn[0];
        printf("structureNumber: %d\n", structureNumber);
        bool useNameAsTitle = (*useNameAsTitleInt!=0);
        fileIter = indigoIterateSDFile(*structures_file);
        int i = 0;
        while ((structure = indigoNext(fileIter))) {
            i++;
            if(i==structureNumber) {
                indigoFoldHydrogens(structure);
                printf("filename: %s", *filename);
                if(useNameAsTitle) {
                    renderMolecule(structure, *filename);
                }
                else {
                    renderMolecule(structure, "", *filename, -1, "");
                }
                break;
            }
            indigoFree(structure);
        }
        indigoFree(fileIter);
    }

    void R_standardiseMolecules(char **structures_file,
                                char **standardised_file,
                                char **removed_file,
                                char **output,
                                int *isSDFInt,
                                int *removeInorganicInt,
                                int *fluorineLimitInt,
                                int *chlorineLimitInt,
                                int *bromineLimitInt,
                                int *iodineLimitInt,
                                int *minMassLimitInt,
                                int *maxMassLimitInt,
                                int *numberProcessedInt) {
        indigoSetOption("aromaticity-model", "generic"); // enable the dearomatization to work properly
        indigoSetOption("dearomatize-verification", "false"); // enable the dearomatization to work properly
        
        bool isSDF = (*isSDFInt!=0);
        bool removeInorganic = (*removeInorganicInt!=0);
        int numberProcessed = *numberProcessedInt;
        int fluorineLimit = *fluorineLimitInt;
        int chlorineLimit = *chlorineLimitInt;
        int bromineLimit = *bromineLimitInt;
        int iodineLimit = *iodineLimitInt;
        int minMassLimit = *minMassLimitInt;
        int maxMassLimit = *maxMassLimitInt;
        
        bool debug = false;
        
        int structure, structureIter;
        int sdfWriter = indigoWriteFile(*standardised_file);
        int removedWriter = indigoWriteFile(*removed_file);
        
        ofstream property_stream;
        
        property_stream.open(*output, ios::out | ios::trunc); // delete the current file
        property_stream.close();
        property_stream.open(*output, ios::out | ios::ate | ios::app | ios::binary);
        
        
        //ofstream property_stream2;
        
        //property_stream2.open("properties.csv", ios::out | ios::trunc); // delete the current file
        //property_stream2.close();
        //property_stream2.open("properties.csv", ios::out | ios::ate | ios::app | ios::binary);

        
        if(isSDF) {
            Rprintf("Reading SDF (C)\n");
            structureIter = indigoIterateSDFile(*structures_file);
        }
        else {
            Rprintf("Reading SMILES (C)\n");
            structureIter = indigoIterateSmilesFile(*structures_file);
        }
        
        int readCount = 0;
        int writeCount = 0;
        int inorganicCount = 0;
        int tooLightCount = 0;
        int tooHeavyCount = 0;
        int highFlourineCount = 0;
        int highChlorineCount = 0;
        int highBromineCount = 0;
        int highIodineCount = 0;
        
        int props,prop,propFirst,propsFirst;
        int propss = 0;
        int namesprop = TRUE;
        property_stream << "NAME" << "\t";
        property_stream << "Kept" << "\n";
        /* property stream contains all properties in the sdf files, whereas property_stream2 only the names and which are kept */
        while (structure = indigoNext(structureIter)) {
            readCount++;
            if(numberProcessed != -1 && readCount > numberProcessed) break;
            // Check if the input structure is Ok.
            if (indigoCountAtoms(structure) == -1) {
                property_stream << readCount << "\t";
                property_stream << "0" << "\n";
            } 
            if (indigoCountAtoms(structure) != -1) {
                int structureIndex = indigoIndex(structure);
                double target_value = 0;
                
                //if(namesprop){
                //    property_stream << "NAME" << "\t"; 
                //    //property_stream2 << "NAME" << "\t";
                //    propsFirst = indigoIterateProperties(structure);
                //    while (propFirst = indigoNext(propsFirst)) {
                //        string prop_val_first = indigoName(propFirst);
                //        property_stream << prop_val_first << "\t";
                //        indigoFree(propFirst);
                //    }
                //    property_stream << "Kept" << "\n";
                //    //property_stream2 << "Kept" << "\n";
                //    namesprop=FALSE;
                //}
                
                
                props = indigoIterateProperties(structure);
                string comp_name = indigoName(structure);
                if (comp_name.size() == 0){ //comp_name.length()
                    std::ostringstream ss;
                    ss << readCount;
                    comp_name = ss.str();
                }
                
                property_stream << comp_name << "\t";
                //property_stream2 << comp_name << "\t";
                //while (prop = indigoNext(props)) {
                //    string target_value = indigoGetProperty(structure, indigoName(prop));
                //    while ( target_value.find ("\n") != string::npos )
                //    {
                //        boost::replace_all(target_value, "\n", ";");
                //    }
                //    property_stream << target_value << "\t";
                //    indigoFree(prop);
                //}
                
                string structureName = indigoName(structure);
                int structureClone = indigoClone(structure);
                
                
                if (debug) Rprintf("%s (#%d)\n", structureName.c_str(), structureIndex+1);
                if (debug) Rprintf("folding hydrogens\n");
                indigoFoldHydrogens(structure);
                
                if((structureIndex+1)%50 == 0) printf(".");
                if((structureIndex+1)%5000 == 0) printf("5000+\n");
                
                if (debug) Rprintf("checking bad valence and ambiguousH\n");
                // skip over if indigo determines bad valence
                if( indigoCheckBadValence(structure)==NULL ) {
                    Rprintf("%s (#%d) skipped over: INDIGO_BAD_VALANCE\n", structureName.c_str(), structureIndex+1);
                    property_stream << "0" << endl;
                    //property_stream2 << "0" << endl;
                    indigoFree(structureClone);
                    indigoFree(structure);
                    continue;
                }
                // skip over if indigo determines ambiguous H
                if( indigoCheckAmbiguousH(structure)==NULL ) {
                    Rprintf("%s (#%d) skipped over: INDIGO_AMBIGUOUSH\n", structureName.c_str(), structureIndex+1);
                    property_stream << "0" << endl;
                    //property_stream2 << "0" << endl;
                    indigoFree(structureClone);
                    indigoFree(structure);
                    continue;
                }
                
                if (debug) Rprintf("checking if organic\n");
                // skip over if not organic
                if(!isOrganic(structure)) {
                    inorganicCount++;
                    if(removeInorganic) {
                        Rprintf("%s (#%d) skipped over: NOT_ORGANIC\n", structureName.c_str(), structureIndex+1);
                        property_stream << "0" << endl;
                        //property_stream2 << "0" << endl;
                        indigoFree(structureClone);
                        indigoFree(structure);
                        continue;
                    }
                }
                
                bool wasAromatic = isAromatic(structure);
                // check if the structure can be dearomatized
                if (debug) Rprintf("dearomitizing\n");
                indigoDearomatize(structure);
                
                if (debug) Rprintf("checking aromatic\n");
                // if the structure is now not aromatic (i.e. hasn't failed the dearomitisation) then pass it through the inchi plugin
                if(!isAromatic(structure)) {
                    if (debug) Rprintf("passing through inchi\n");
                    const char* inchi = indigoInchiGetInchi(structure);
                    string trimmedInchi = trimInchi(inchi);
                    int temp = indigoInchiLoadMolecule(trimmedInchi.c_str());
                    if(strcmp(indigoInchiGetWarning(), "") != 0) {
                        Rprintf("%s (#%d) warning: while converting to Inchi: %s\n", structureName.c_str(), structureIndex+1, indigoInchiGetWarning());

                    }

                    else {
                       indigoFree(structure);
                       structure = temp;
                    }
                }
                else {
                    if(wasAromatic) {
                        Rprintf("%s (#%d) warning: failed dearomatize, didn't go through inchi\n", structureName.c_str(), structureIndex+1);

                    }
                }
                
                
                // reset the isotopes of all the atoms in the molecule
                if (debug) Rprintf("resetting isotopes\n");
                resetIsotopes(structure);
                
                // print warnings on molecular mass variations
                if (debug) Rprintf("various conditions\n");
                float molecularMass = indigoMolecularWeight(structure);
                if(minMassLimit != -1 && molecularMass < minMassLimit) {
                    tooLightCount++;
                    Rprintf("%s (#%d) warning: molecular mass less than %d daltons\n", structureName.c_str(), structureIndex+1, minMassLimit);
                }
                if(maxMassLimit != -1 && molecularMass > maxMassLimit) {
                    tooHeavyCount++;
                    Rprintf("%s (#%d) warning: molecular mass greater than %d daltons\n", structureName.c_str(), structureIndex+1, maxMassLimit);
                }
                bool highFlourine = fluorineLimit != -1 && containsMoreThanX(structure, fluorineLimit, 9);
                if(highFlourine) {
                    highFlourineCount++;
                    Rprintf("%s (#%d) warning: molecule contains more than %d flourines\n", structureName.c_str(), structureIndex+1, fluorineLimit);
                }
                bool highChlorine = chlorineLimit != -1 && containsMoreThanX(structure, chlorineLimit, 17);
                if(highChlorine) {
                    highChlorineCount++;
                    Rprintf("%s (#%d) warning: molecule contains more than %d chlorines\n", structureName.c_str(), structureIndex+1, chlorineLimit);
                }
                bool highBromine = bromineLimit != -1 && containsMoreThanX(structure, bromineLimit, 35);
                if(highBromine) {
                    highBromineCount++;
                    Rprintf("%s (#%d) warning: molecule contains more than %d bromines\n", structureName.c_str(), structureIndex+1, bromineLimit);
                }
                bool highIodine = iodineLimit != -1 && containsMoreThanX(structure, iodineLimit, 53);
                if(highIodine) {
                    highIodine++;
                    Rprintf("%s (#%d) warning: molecule contains more than %d iodines\n", structureName.c_str(), structureIndex+1, iodineLimit);
                }
                
                // remove 'bad' molecules (normally used during training)
                if(maxMassLimit != -1){
                    if(molecularMass < minMassLimit || molecularMass > maxMassLimit || highFlourine || highChlorine || highBromine || highIodine) {
                        indigoSdfAppend(removedWriter, structure);
                        indigoFree(structureClone);
                        indigoFree(structure);
                        property_stream << "0" << endl;
                        continue;
                    } else {
                        
                        property_stream << "1" << endl;
                        
                    }
                } else {
                    
                    if(molecularMass < minMassLimit || highFlourine || highChlorine || highBromine || highIodine) {
                        indigoSdfAppend(removedWriter, structure);
                        indigoFree(structureClone);
                        indigoFree(structure);
                        property_stream << "0" << endl;
                        continue;
                    } else {
                        
                        property_stream << "1" << endl; 
                    }
                }
                
                
                // use the largest substructure
                if (debug) Rprintf("selecting the largest substructure\n");
                int temp2 = pickLargestSubstructure(structure);
                indigoFree(structure);
                structure = temp2;
                
                writeCount++;
                std::string addition = "Standardised_";
                std::string Mname = indigoName(structureClone);
                if (Mname.size() == 0){
                    std::ostringstream ss;
                    ss << readCount;
                    Mname = ss.str();
                }
                
                indigoSetName(structure, (addition + Mname).c_str());
///
//add fields to the sdf file
               props = indigoIterateProperties(structureClone);
               while (prop = indigoNext(props)) {
                   string prop_name = indigoName(prop);
				   string prop_value = indigoGetProperty(structureClone, indigoName(prop));
				   indigoSetProperty(structure, prop_name.c_str(), prop_value.c_str());
                   indigoFree(prop);
                    }

////

                indigoSdfAppend(sdfWriter, structure);
                indigoFree(structureClone);
                indigoFree(structure);
            }
            else {
                Rprintf("%s, Readcount: %d\n", indigoGetLastError(), readCount);
            }
        }
        printf("\n");
        indigoFree(structureIter);
        indigoClose(sdfWriter);
        indigoFree(sdfWriter);
        indigoClose(removedWriter);
        indigoFree(removedWriter);
        property_stream.close();
        
        printf("\ninorganicCount: %d\n", inorganicCount);
        printf("tooLightCount: %d\n", tooLightCount);
        printf("tooHeavyCount: %d\n", tooHeavyCount);
        printf("highFlourineCount: %d\n", highFlourineCount);
        printf("highChlorineCount: %d\n", highChlorineCount);
        printf("highBromineCount: %d\n", highBromineCount);
        printf("highIodineCount: %d\n", highIodineCount);
        printf("readCount: %d\n", readCount);
        printf("writeCount: %d\n", writeCount);
    }

} // extern "C"
    
    
    
    
    
    
    
