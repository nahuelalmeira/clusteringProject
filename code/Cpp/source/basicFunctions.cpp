#include "../headers/basicFunctions.hpp"
#include <string>

using std::string;

void createDir(string dir_name) {
    string command = "mkdir -p " + dir_name;
    system(command.c_str());
}

/* Compute some exponential values in advance for increasing 
   performance in the MC simulation 
   TODO: not implemented yet */
/*
void computeExponentials(double mu) {
    if(mu<0) mu = -mu; //TODO: test that this is ok for negative mu
	for(int i=0; i < EXP_VALUES; i++)
		exponential[i] = exp(-mu*i);
}
*/