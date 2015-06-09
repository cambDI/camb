void renderMolecule(int structure, const char* folder, const char* file, int number, const char* title);

void renderMolecule(int structure, const char* file);

void renderPair(int s1, int s2, const char* folder, const char* file, int number, const char* title);

int neutraliseCSPnegO(int structure);

int neutraliseNposH(int structure);

int pickLargestSubstructure(int structure);

bool isOrganic(int structure);

bool isAromatic(int structure);

void resetIsotopes(int structure);

// returns true if the structure has more than X atoms with atomic number Y
bool containsMoreThanX(int structure, int X, int atomicNumber);
