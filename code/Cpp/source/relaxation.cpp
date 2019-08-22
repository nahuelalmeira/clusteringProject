/*

Parameters:

--randMCS: number of Monte Carlo steps to randomize the network
--connected: if present, rewirings will take place only if they don't disconnect the network
--seed: random seed
--samples: number of samples 
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
#include <vector>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::setprecision;
using std::string;
using std::to_string;
using std::fstream;
using std::ios;
using std::swap;
using std::vector;

string NETWORK;
string MODE = "random";
long seed = -17379;       // Arbitrary negative number for ran2
string str_seed = "0";
int RAND_MCS = 100; 
int SAMPLES = 1000;     
bool CONNECTED = false;
bool VERBOSE = false;
   
string INPUT_DIR = "../networks";
string INPUT_FILE;

vector<double> get_log_values(int samples, int rand_mcs) {
    vector<double> log_values;

    double min_val = 0.01;
    double max_val = (double) rand_mcs;
    double log_min_val = log10(min_val);
    double log_max_val = log10(max_val);
    double step = (log_max_val - log_min_val) / samples;
    for(int i = 0; i < samples; i++)
        log_values.push_back(pow(10., log_min_val + i*step));
    for(int i = 0; i < samples; i++)
        log_values.push_back((double) (samples + i));
    return log_values;
}

/* Parse command line arguments */
void parseArguments(int argc, char *argv[]) {

    NETWORK = argv[1];
    for(int i=2; i<argc; i++) {
        if ( !strcmp(argv[i], "--seed") ) {
            seed += atoi(argv[i+1]);
            str_seed = argv[i+1];
        }
        if ( !strcmp(argv[i], "--connected") )
            CONNECTED = true;
        if( !strcmp(argv[i], "--inputDir") )
            INPUT_DIR = argv[i+1];
        if( !strcmp(argv[i], "--inputFile") )
            INPUT_FILE = argv[i+1];
        if ( !strcmp(argv[i], "--randMCS") )
            RAND_MCS = atoi(argv[i+1]);
        if ( !strcmp(argv[i], "--samples") )
            SAMPLES = atoi(argv[i+1]);
        if ( !strcmp(argv[i], "--verbose") )
            VERBOSE = true;
    }
    if(INPUT_FILE == "") INPUT_FILE = NETWORK + ".txt";
}

/* Print arguments */
void printArguments(ostream& out) {
    out << "\n---- Arguments ----" << endl;
    out << "Network     = " << NETWORK << endl;
    out << "Input file  = " << INPUT_FILE << endl;
    out << "mode        = " << MODE << endl;
    out << "randMCS     = " << RAND_MCS << endl;
    out << "seed        = " << seed << " + " << str_seed << endl;
    out << endl;
}

string protocolName() {
    string str_randMCS = to_string(RAND_MCS);
    string str_samples = to_string(SAMPLES);
    string protocol;
    if (CONNECTED) protocol = "relaxation_connected_randMCS" + str_randMCS + "_samples" + str_samples;
    else           protocol = "relaxation_randMCS" + str_randMCS + "_samples" + str_samples;
    return protocol;
}

string buildOutputDir() {
    string output_dir;
    string protocol = protocolName();

    output_dir = INPUT_DIR + "/" + NETWORK;
    //if(CONNECTED) output_dir += "/connected";
    output_dir += "/" + protocol + "/seed" + string(5 - str_seed.length(), '0') + str_seed; 
    return output_dir;
}

int main(int argc, char *argv[]) {

    int valid_swap = 0;
    long valid_swaps = 0;
    double C, Cws;
    double delta_C = 0;
    double delta_Cws = 0;
    Edge e1, e2;

    parseArguments(argc, argv);
    printArguments(cout);

    vector<double> log_values = get_log_values(SAMPLES, RAND_MCS);

    string edgelist_file_name;

    /* Build graphs */
    Graph G_original = createGraph(INPUT_DIR + "/" + NETWORK + "/" + INPUT_FILE, true);
    Graph G = createGraph(INPUT_DIR + "/" + NETWORK + "/" + INPUT_FILE, false);

    G.summary(cout);

    C = G.C;
    Cws = G.Cws;

    string output_dir = buildOutputDir();
    createDir(output_dir);

    valid_swaps = 0;
    int k = 0;
    int checkpoint = 0;
    double log_value;
    string str_log_value;
    //cout << "log_values = " << endl;
    //for(int i = 0; i<log_values.size(); i++)
    //    cout << log_values[i] << endl;
    
    //while(valid_swaps < RAND_MCS*G.M) {
    while(k < int(log_values.size())) {

        // Perform a valid MCS
        do {
            selectPairOfLinks(G, &e1, &e2, &seed);
            valid_swap = doubleEdgeSwap(G, e1, e2, &delta_C, &delta_Cws, 
                                        "random", CONNECTED, &seed, 0);
        } while(!valid_swap);
        ++valid_swaps;

        // Update clustering values
        C += delta_C;
        Cws += delta_Cws;
        G.C = C;
        G.Cws = Cws;

        // Write data to edgelist and print to standard output
        if(valid_swaps >= checkpoint) {
            //str_sample = string( - to_string(M).length(), '0') + to_string(M);
            log_value = log_values[k++];
            str_log_value = to_string(log_value);
            str_log_value = string(11 - to_string(log_value).length(), '0') + to_string(log_value);
            edgelist_file_name = output_dir + "/" "M" + str_log_value + ".txt";
            G.writeEdgeList(edgelist_file_name);
            checkpoint = log_value*G.M;
            if(VERBOSE) cout << k << " " << log_value << " " << G.C << " " << G.Cws << endl;
        }          
    }

    G.summary(cout);
    testConsistency(G);

    cout << "original G connected = " << G_original.isConnected() << endl;
    cout << "rewired  G connected = " << G.isConnected() << endl;
    
    return 0;
}