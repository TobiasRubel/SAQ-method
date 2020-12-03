#if ! defined (TYPEDEFS_H) 
#define TYPEDEFS_H

/* the headers and typedefs used throughout the code */

#include <map>
#include <vector>
#include <iostream>
#include <armadillo>

//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#endif

using namespace std;
using namespace arma;

typedef map<string,string> SShash;
typedef map<string,int> SIhash;
typedef vector<int> split;

typedef vector<double> VF;
typedef vector<VF> MF;
typedef vector<MF> T3F;
typedef vector<T3F> T4F;
typedef vector<T4F> vec_tensor;
typedef vector<mat> VM;
typedef vector<VM> MM;

extern int VERBOSITY;
extern unsigned int NUM_SPECIES;

#endif
