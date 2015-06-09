ifeq ($(shell uname), Linux)
	INDIGO_LIBRARIES := -lfontconfig -lfreetype
endif

all:clean
	PKG_CXXFLAGS="-I$(PWD)/api -I$(PWD)/api/plugins/renderer -I$(PWD)/api/plugins/inchi"  R CMD SHLIB -o camb.so GetProperties.cpp fingerprint.cpp standardise.cpp misc.cpp smindigo.cpp -I$(PWD)/api -I$(PWD)/api/plugins/renderer -I$(PWD)/api/plugins/inchi -L$(PWD)/build/dist/Mac/10.6/lib -L$(PWD)/build/dist/Mac/10.7/lib -L$(PWD)/build/dist/Mac/10.8/lib -L$(PWD)/build/dist/Mac/10.9/lib -L$(PWD)/build/dist/Mac/10.10/lib -L$(PWD)/build/dist/Linux/x64/lib -L$(PWD)/build/dist/Linux/x86/lib  -L$(PWD)/build/dist/Win/x64/lib -L$(PWD)/build/dist/Win/x86/lib -lindigo-static -lcommon -lmolecule -lreaction -llayout -lz -ltinyxml -lgraph -lindigo-inchi-static -lindigo-renderer-static -lpng -lpixman -lcairo -linchi -lrender2d -lindigo-static -lcommon -lmolecule -lreaction -llayout -lz -ltinyxml -lgraph -lindigo-inchi-static -lindigo-renderer-static -lpng -lpixman -lcairo -linchi -lrender2d -lindigo-static -lcommon -lmolecule -lreaction -llayout -lz -ltinyxml -lgraph -lindigo-inchi-static -lindigo-renderer-static -lpng -lpixman -lcairo -linchi -lrender2d $(INDIGO_LIBRARIES)
clean:
	rm -f *.o
	rm -f *.so

