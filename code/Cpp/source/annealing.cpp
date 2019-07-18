/*

Parameters:

--randMCS: number of Monte Carlo steps to randomize the network
--connected: if present, rewirings will take place only if they don't disconnect the network
--mode: maxC, maxCws
--seed: random seed
--initialMu:
--finalMu:
--deltaMu:
--cycle: change deltaMu by -deltaMu and continue until mu == initialMu
--samples: number of samples for each mu
--decorrTime: time between samples
--transitTime: time before first sample for each mu
--inputDir: input directory
--inputFile: input file, relative to inputDir
*/

#include "../headers/basicFunctions.hpp"
#include "../headers/graphFunctions.hpp"
#include <iostream>
#include <string>
#include <string.h>
#include <iomanip>           // setprecision
#include <sstream>           // stringstream

using std::cout;
using std::endl;
using std::ostream;
using std::setprecision;
using std::string;
using std::to_string;
using std::fstream;
using std::ios;
using std::swap;

//double exponential[EXP_VALUES];

string NETWORK;
string MODE = "random";
long seed = -17379;       // Arbitrary negative number for ran2
string str_seed = "0";
int RAND_MCS = 100;     
bool CONNECTED = false;
double INITIAL_MU = 0.;
double FINAL_MU = 10.;
double DELTA_MU = 0.1; 
bool CYCLE = false;   
string INPUT_DIR = "../networks";
string INPUT_FILE;
string RAMPE = "annealing";
int N_SAMPLES = 100;
int DECORR_TIME = 100;
int TRANSIT_TIME = 1000;


/* Parse command line arguments */
void parseArguments(int argc, char *argv[]) {

    NETWORK = argv[1];
    for(int i=2; i<argc; i++) {
        if ( !strcmp(argv[i], "--mode") || !strcmp(argv[i], "-m") )
            MODE = argv[i+1];
        if ( !strcmp(argv[i], "--seed") ) {
            seed += atoi(argv[i+1]);
            str_seed = argv[i+1];
        }
        if ( !strcmp(argv[i], "--randMCS") )
            RAND_MCS = atoi(argv[i+1]);
        if ( !strcmp(argv[i], "--connected") )
            CONNECTED = true;
        if( !strcmp(argv[i], "--initialMu") )
            INITIAL_MU = atof(argv[i+1]);
        if( !strcmp(argv[i], "--finalMu") )
            FINAL_MU = atof(argv[i+1]);
        if( !strcmp(argv[i], "--deltaMu") ) {
            DELTA_MU = atof(argv[i+1]);
            if(DELTA_MU > 0) RAMPE = "annealing";
            else             RAMPE = "cooling";
        }
        if( !strcmp(argv[i], "--cycle") )
            CYCLE = true;
        if( !strcmp(argv[i], "--inputDir") )
            INPUT_DIR = argv[i+1];
        if( !strcmp(argv[i], "--inputFile") )
            INPUT_FILE = argv[i+1];

        if( !strcmp(argv[i], "--decorrTime") )
            DECORR_TIME = atoi(argv[i+1]);
        if( !strcmp(argv[i], "--transitTime") )
            TRANSIT_TIME = atoi(argv[i+1]);
        if( !strcmp(argv[i], "--samples") )
            N_SAMPLES = atoi(argv[i+1]);
    }
    if(INPUT_FILE == "") INPUT_FILE = NETWORK + ".txt";
    //if( !RAMPE.compare("cooling") )
    if (fabs(INITIAL_MU) > 10E-6)
        RAND_MCS = 0;

}

/* Print arguments */
void printArguments(ostream& out) {
    out << "\n---- Arguments ----" << endl;
    out << "Network     = " << NETWORK << endl;
    out << "Input file  = " << INPUT_FILE << endl;
    out << "mode        = " << MODE << endl;
    out << "randMCS     = " << RAND_MCS << endl;
    out << "seed        = " << seed << " + " << str_seed << endl;
    out << "initialMu   = " << INITIAL_MU << endl;
    out << "finalMu     = " << FINAL_MU << endl;
    out << "deltaMu     = " << DELTA_MU << endl;
    out << "cycle       = " << CYCLE << endl;
    out << "samples     = " << N_SAMPLES << endl;
    out << "transitTime = " << TRANSIT_TIME << endl;
    out << "decorrTime  = " << DECORR_TIME << endl; 
    out << endl;
}

/* Write data to file */
void writeData(const Graph G, double mu, int sample, int steps, ostream& out) {

    out << mu  << " " << sample << " "         
        << G.C << " " << G.Cws  << " "
        << steps  << endl;
}

string protocolName(string mode, string rampe, double mu_init, double step, int samples, int transit, int decorr) {
    string protocol;

    string str_decorr;
    string str_transit;
    string str_samples;
    string str_step;  
    string str_mu_init;

    str_transit = string(6 - to_string(transit).length(), '0') + to_string(transit);
    str_decorr  = string(6 - to_string(decorr).length(), '0')  + to_string(decorr);
    str_samples = string(6 - to_string(samples).length(), '0') + to_string(samples);
    str_mu_init = to_string(mu_init);
    str_step = to_string(fabs(step));

    protocol = mode + "_" + rampe + "_muInit" + str_mu_init + "_step" + str_step 
             + "_samples" + str_samples + "_transit" + str_transit + "_decorr" + str_decorr;

    return protocol;
}

string buildOutputDir(string mode, string rampe, double mu_init, double step, int samples, int transit, int decorr) {
    string output_dir;
    string protocol = protocolName(mode, rampe, mu_init, step, samples, transit, decorr);

    output_dir = INPUT_DIR + "/" + NETWORK;
    if(CONNECTED) output_dir += "/connected";
    output_dir += "/" + protocol + "/seed" + string(5 - str_seed.length(), '0') + str_seed; 
    return output_dir;
}

int main(int argc, char *argv[]) {

    int valid_swap = 0;
    long valid_swaps = 0;
    long total_it = 0;
    double C, Cws;
    double delta_C = 0;
    double delta_Cws = 0;
    Edge e1, e2;
    int sample = 1;

    parseArguments(argc, argv);
    printArguments(cout);

    string edgelist_file_name;

    /* Build graphs */
    Graph G_original = createGraph(INPUT_DIR + "/" + NETWORK + "/" + INPUT_FILE, true);
    Graph G = createGraph(INPUT_DIR + "/" + NETWORK + "/" + INPUT_FILE, false);

    G.summary(cout);

    C = G.C;
    Cws = G.Cws;

    /* Randomize network by double swapping edges */
    if (RAND_MCS) cout << "Randomizing network with mu = 0." << endl;
    while(valid_swaps < RAND_MCS*G.M) {
        
        // Perform a valid MCS
        do {
            selectPairOfLinks(G, &e1, &e2, &seed);
            valid_swap = doubleEdgeSwap(G, e1, e2, &delta_C, &delta_Cws, 
                                        "random", CONNECTED, &seed, 0);
            total_it++;
        } while(!valid_swap);
        ++valid_swaps;

        // Update clustering values
        C += delta_C;
        Cws += delta_Cws;
        G.C = C;
        G.Cws = Cws;
    }
    cout << "End randomization." << endl;

    cout << "original G connected = " << G_original.isConnected() << endl;
    cout << "rewired  G connected = " << G.isConnected() << endl;

    total_it = 0;
    int total_steps;
    int steps;
    bool write = false;
    double mu = INITIAL_MU;

    string str_sample;
    string str_steps;

    string output_dir = buildOutputDir(MODE, RAMPE, INITIAL_MU, DELTA_MU, N_SAMPLES, TRANSIT_TIME, DECORR_TIME);
    createDir(output_dir);

    while(true) {

        cout << "mu = " << setprecision(3) << mu << endl;

        sample = 1;
        steps = 0;
        valid_swaps = 0;
        total_steps = TRANSIT_TIME + (N_SAMPLES-1)*DECORR_TIME;
        while(steps <= total_steps) {

            // Perform a valid MCS
            do {
                selectPairOfLinks(G, &e1, &e2, &seed);
                valid_swap = doubleEdgeSwap(G, e1, e2, &delta_C, &delta_Cws, 
                                            MODE, CONNECTED, &seed, mu);
            } while(!valid_swap);
            ++valid_swaps;

            // Update clustering values
            C += delta_C;
            Cws += delta_Cws;
            G.C = C;
            G.Cws = Cws;

            // Write data to edgelist and print to standard output
            if(valid_swaps == G.M) {
                valid_swaps = 0;
                steps++;
                write = true;
            }
            if(steps >= TRANSIT_TIME) {
                if(steps%DECORR_TIME==0 && write)  {
                    str_sample = string(4 - to_string(sample).length(), '0') + to_string(sample);
                    str_steps = string(6 - to_string(steps).length(), '0') + to_string(steps);
                    edgelist_file_name = output_dir + "/" + 
                                         RAMPE + "_mu" + to_string(mu) + "_sample" + str_sample + ".txt";
                    G.writeEdgeList(edgelist_file_name);
                    write = false;
                    sample++; 
                }  
            }
        }

        // Update mu
        mu += DELTA_MU;
        if( fabs(mu-FINAL_MU) < 10E-6 ) {
            if(!CYCLE) break;
            else { // Now the reverse simulation will start
                swap(INITIAL_MU, FINAL_MU);
                DELTA_MU = -DELTA_MU;
                CYCLE = false;
                if      (RAMPE == "annealing") RAMPE = "cooling";
                else if (RAMPE == "cooling")   RAMPE = "annealing";
                output_dir = buildOutputDir(MODE, RAMPE, INITIAL_MU, DELTA_MU, N_SAMPLES, TRANSIT_TIME, DECORR_TIME);
                createDir(output_dir);
            }
        }
    }

    G.summary(cout);
    testConsistency(G);

    cout << "original G connected = " << G_original.isConnected() << endl;
    cout << "rewired  G connected = " << G.isConnected() << endl;
    
    return 0;
}