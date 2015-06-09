#include <stdio.h>
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

//#include <string.hpp>
//#include <boost/algorithm/string.hpp>

#include <indigo.h>
#include "indigo-renderer.h"
#include "indigo-inchi.h"
#include "smindigo.h"

using namespace std;


/// todo
///rm properties.csv file with sys
/// create vector when reading only one property
/// check for smiles

//// plotting molecules


extern "C" {
void ReplaceStringInPlace(std::string& subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
          pos += replace.length();
      }
}





    /// Get all properties in structures file
    void R_GetPropertiesSDF(char **structures_file,
                            //  char **standardised_file,
                            //char **target_field_name,
                            int *numberProcessedInt,
           //                 int *isSDFInt,
                            char **output) {
        indigoSetOption("aromaticity-model", "generic"); // enable the dearomatization to work properly
        indigoSetOption("dearomatize-verification", "false"); // enable the dearomatization to work properly
        //bool isSDF = (*isSDFInt!=0);
        int numberProcessed = *numberProcessedInt;
        
        bool debug = false;
        
        int structure, structureIter,props,prop,propFirst,propsFirst;
        ///int sdfWriter = indigoWriteFile(*standardised_file);
        /// int removedWriter = indigoWriteFile(*removed_file);
        
        ofstream target_stream;
        
        Rprintf("Properties will be written");
        target_stream.open(*output, ios::out | ios::trunc); // delete the current file
        target_stream.close();
        target_stream.open(*output, ios::out | ios::ate | ios::app | ios::binary);
        
        
        //if(isSDF) {
            Rprintf("Reading SDF (C)\n");
            structureIter = indigoIterateSDFile(*structures_file);
        //}
        //else {
          //  Rprintf("Reading SMILES (C)\n");
          //  structureIter = indigoIterateSmilesFile(*structures_file);
       // }
        
        int readCount = 0;
        int writeCount = 0;
        
        int propss = 0;
        int namesprop = 0;
        
        while (structure = indigoNext(structureIter)) {
            readCount++;
            if(numberProcessed != -1 && readCount > numberProcessed) break;
            if (indigoCountAtoms(structure) != -1) {
                int structureIndex = indigoIndex(structure);
                double target_value = 0;
                //if(write_targets) {
                    
                    /// we write property names
                    
                    if(namesprop == 0){
                        propsFirst = indigoIterateProperties(structure);
                        while (propFirst = indigoNext(propsFirst)) {
                            string prop_val_first = indigoName(propFirst);
                            while ( prop_val_first.find ("\n") != string::npos )
                            {
                                //boost::replace_all(prop_val_first, "\n", ";");
                                ReplaceStringInPlace(prop_val_first, ";", "\n");
                            }
                            
                            target_stream << prop_val_first << "\t";
                            indigoFree(propFirst);
                        }
                        target_stream << "\n";
                        namesprop=1;
                    }
                    
                    
                    props = indigoIterateProperties(structure);
                    
                    while (prop = indigoNext(props)) {
                        string prop_value = indigoGetProperty(structure, indigoName(prop));
                        while ( prop_value.find ("\n") != string::npos )
                        {
                            //boost::replace_all(prop_value, "\n", ";");
                            ReplaceStringInPlace(prop_value, ";", "\n");
                        }
                        target_stream << prop_value << "\t";
                        indigoFree(prop);
                    }
                    target_stream << "\n";

                string structureName = indigoName(structure);
                printf("Structure Index: %d\n", structureIndex+1);
                
                writeCount++;
                
                indigoFree(structure);
            }
            else {
                Rprintf("%s, Readcount: %d\n", indigoGetLastError(), readCount);
            }
        }
        printf("\n");
        indigoFree(structureIter);
        target_stream.close();
        
        
        printf("readCount: %d\n", readCount);
    }
    
    
    
    
    /// Get unique property
    
    void R_GetPropertySDF(char **structures_file,
                          char **target_field_name,
                          int *numberProcessedInt,
                          //int *isSDFInt,
                          char **output) {
        indigoSetOption("aromaticity-model", "generic"); // enable the dearomatization to work properly
        indigoSetOption("dearomatize-verification", "false"); // enable the dearomatization to work properly
        //bool isSDF = (*isSDFInt!=0);
        bool write_targets = !(*target_field_name && *target_field_name[0] == '\0');
        int numberProcessed = *numberProcessedInt;
        
        bool debug = false;
        int structure, structureIter;
        ofstream target_stream;
        
        if(write_targets) {
            Rprintf("Properties will be written");
            target_stream.open(*output, ios::out | ios::trunc); // delete the current file
            target_stream.close();
            target_stream.open(*output, ios::out | ios::ate | ios::app | ios::binary);
        }
        
        //if(isSDF) {
            Rprintf("Reading SDF (C)\n");
            structureIter = indigoIterateSDFile(*structures_file);
        //}
        //else {
            //Rprintf("Reading SMILES (C)\n");
            //structureIter = indigoIterateSmilesFile(*structures_file);
        //}
        
        int readCount = 0;
        int writeCount = 0;
        int namesprop = 0;
        
        while (structure = indigoNext(structureIter)) {
            readCount++;
            if(numberProcessed != -1 && readCount > numberProcessed) break;
            if (indigoCountAtoms(structure) != -1) {
                int structureIndex = indigoIndex(structure);
                double target_value = 0;
                if(write_targets) {

                    if(indigoHasProperty(structure, *target_field_name)) {
                        if(namesprop == 0){
                            string prop_val_first = indigoGetProperty(structure, *target_field_name);
                            while ( prop_val_first.find ("\n") != string::npos )
                            {
                                //boost::replace_all(prop_val_first, "\n", ";");
                                ReplaceStringInPlace(prop_val_first, ";", "\n");
                            }
                            target_stream << prop_val_first << "\n";
                            namesprop=1;
                        }
                        string prop_value = indigoGetProperty(structure, *target_field_name);
                        while ( prop_value.find ("\n") != string::npos )
                        {
                            //boost::replace_all(prop_value, "\n", ";");
                            ReplaceStringInPlace(prop_value, ";", "\n");
                        }
                        target_stream << prop_value << endl;
                    }
                    else {
                        write_targets = false;
                        target_stream << "no property found on first molecule in file with fieldname: " << *target_field_name << endl;
                        target_stream << " Thus, no the calculation will be stopped.";
                        break;
                    }
                    
                }
                printf("Structure Index: %d\n", structureIndex+1);
                writeCount++;
                indigoFree(structure);
            }
            else {
                Rprintf("%s, Readcount: %d\n", indigoGetLastError(), readCount);
            }
        }
        printf("\n");
        indigoFree(structureIter);
        target_stream.close();
        printf("readCount: %d\n", readCount);
    }
    
    
    /// Show properties of molecules file
    
    void R_ShowPropertiesSDF(char **structures_file,
                            char **output//,
                           /// int *numberProcessedInt,
                          ///int *isSDFInt
						  ) {
        indigoSetOption("aromaticity-model", "generic"); // enable the dearomatization to work properly
        indigoSetOption("dearomatize-verification", "false"); // enable the dearomatization to work properly
        //bool isSDF = (*isSDFInt!=0);
       /// int numberProcessed = *numberProcessedInt;
        
        bool debug = false;
        int structure, structureIter,propFirst,propsFirst;
        ofstream target_stream;
        
            target_stream.open(*output, ios::out | ios::trunc); // delete the current file
            target_stream.close();
            target_stream.open(*output, ios::out | ios::ate | ios::app | ios::binary);
        
        //if(isSDF) {
            Rprintf("Reading SDF (C)\n");
            structureIter = indigoIterateSDFile(*structures_file);
        //}
        //else {
//            Rprintf("Reading SMILES (C)\n");
//            structureIter = indigoIterateSmilesFile(*structures_file);
        //}
        
        int namesprop = 0;
        
        while (structure = indigoNext(structureIter)) {
            if (indigoCountAtoms(structure) != -1) {
                double target_value = 0;
                
                if(namesprop == 0){
                    propsFirst = indigoIterateProperties(structure);
                    while (propFirst = indigoNext(propsFirst)) {
                        string prop_val_first = indigoName(propFirst);
                        while ( prop_val_first.find ("\n") != string::npos )
                        {
                            //boost::replace_all(prop_val_first, "\n", ";");
                            ReplaceStringInPlace(prop_val_first, ";", "\n");
                        }
                        target_stream << prop_val_first << endl;
                        indigoFree(propFirst);
                    }
                    ///target_stream << "\n";
                    namesprop=1;
                    break;
                }
                    else {
                        target_stream << "Incorrect molecules.." << endl;
                        target_stream << " Thus, no the calculation will be stopped.";
                        break;
                    }
            
            indigoFree(structure);
            }
	
            else {
                Rprintf("%s\n", indigoGetLastError());
            }
        }
        indigoFree(structureIter);
        target_stream.close();
    }
    
    
    
    
    
    
} // extern "C"







