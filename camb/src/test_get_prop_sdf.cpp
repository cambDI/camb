//
//  test_get_prop_sdf.cpp
//  
//
//  Created by Isidro on 11/24/13.
//
//

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include "indigo.h"
#include "indigo-renderer.h"
#include "indigo-inchi.h"
#include "./smindigo.h"

using namespace std;



extern "C" {
    
///    void R_getPropertySDF(char **structures_file, int *structureNumberIn, char **property, int *useNameAsTitleInt) {
///        int structure, fileIter;
///		int structureNumber = structureNumberIn[0];
///        bool useNameAsTitle = (*useNameAsTitleInt!=0);
///        fileIter = indigoIterateSDFile(*structures_file);
///        int i = 0;
///        while ((structure = indigoNext(fileIter))) {
///            i++;
///            if(i==structureNumber) {
///                if (indigoHasProperty(structure, *property)) {
///                    indigoGetProperty(structure, *property);
///					printf("aa");
///                }
///            indigoFree(structure);
///		    } 
///		}
///        indigoFree(fileIter);
///       
///	}
///


    void R_getPropertySDF(char **structures_file, char **property, int *useNameAsTitleInt) {
        int structure, fileIter;
        bool useNameAsTitle = (*useNameAsTitleInt!=0);
        fileIter = indigoIterateSDFile(*structures_file);
        int i = 0;
        while ((structure = indigoNext(fileIter))) {
            i++;
                if (indigoHasProperty(structure, *property)) {
                    indigoGetProperty(structure, *property);
					printf("aa");
                }
            indigoFree(structure);
		}
        indigoFree(fileIter);
       
	}






} // extern "C"
