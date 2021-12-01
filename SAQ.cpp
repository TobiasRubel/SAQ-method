// SAQ.cpp V.1.0
// Date: Dec/2020

// This file is part of the "SAQ" program
/*
This code has been written by  M. Garrote-L�pez. 
Please quote the paper: M. Casanellas, J. Fern�ndez-S�nchez, M. Garrote-L�pez
"SAQ: semi-algebraic quartet reconstruction method", https://arxiv.org/abs/2011.13968
*/

#include "functionsSAQ.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;
using namespace arma; // for 'armadillo';

const string s[] = {"A", "C", "G", "T"};
const double epsilon = 1e-16;

int main(int argc, char **argv)
{
    // INPUT PATH
    if (argc < 2)
    {
        cerr << "Missing filename" << endl;
        return 1;
    }
    string filename = argv[1];

    string outfilename;
    if (argv[2] == NULL)
    {
        outfilename = "output.txt";
    }
    else
    {
        outfilename = argv[2];
    }

    double filter;
    if (argc < 3 || argv[3] == NULL)
    {
        filter = -1;
    }
    else
    {
        filter = stod(argv[3]);
    }

    double num_threads;
    if (argc < 4 || argv[4] == NULL) 
    {
        num_threads = 1;
    }
    else
    {
        num_threads = stod(argv[4]);
    }
    omp_set_num_threads(num_threads);
    cout << "Running SAQ with " << num_threads << " threads.\n";

    // READ DATA
    Alignment align = readFASTA(filename);
    int length = align.seq_len;

    vector<Alignment> subsets = GetAlignmentSubsets(align);
    cout << "Number of subsets: " << subsets.size() << "\n";

    string bestQuartets[subsets.size()];

    #pragma omp parallel 
    {
    #pragma omp for
    for (int i = 0; i < subsets.size(); ++i) {
        map<string, int> tensor = Getcolumns(subsets[i]);

        // SAQ method
        T4F Mtensor = convert_tensor_4(tensor, length);
        T4F zeroTens = T4F(4, T3F(4, MF(4, VF(4, 0))));
        MM Ns = double_marginalizations(Mtensor);
        int rank_matrix = 4;
        double eps = 1e-7;

        vec aux_01 = SAQ(Mtensor, {0, 1, 2, 3}, Ns, rank_matrix, filter);
        vec aux_02 = SAQ(Mtensor, {0, 2, 1, 3}, Ns, rank_matrix, filter);
        vec aux_03 = SAQ(Mtensor, {0, 3, 1, 2}, Ns, rank_matrix, filter);

        double score_01, score_02, score_03;
        double weight_01, weight_02, weight_03;

        score_01 = aux_01[0];
        score_02 = aux_02[0];
        score_03 = aux_03[0];

        double sum = score_01 + score_02 + score_03;
        string sp0, sp1, sp2, sp3;
        sp0 = subsets[i].taxa[0];
        sp1 = subsets[i].taxa[1];
        sp2 = subsets[i].taxa[2];
        sp3 = subsets[i].taxa[3];
        
        if (score_02 >= score_01 && score_02 >= score_03)
        {
            bestQuartets[i] = sp0 + "," + sp2 + "|" + sp1 + "," + sp3;
        }
        else if (score_03 >= score_01 && score_03 >= score_01)
        {
            bestQuartets[i] = sp0 + "," + sp3 + "|" + sp1 + "," + sp2;
        }
        else //if (score_01 >= score_02 && score_01 >= score_03)
        {
            bestQuartets[i] = sp0 + "," + sp1 + "|" + sp2 + "," + sp3;
        }
    }
    }

    string quartetOutput = "";
    for (int i = 0; i < subsets.size(); i++) {
        quartetOutput += bestQuartets[i] + " ";
    }

    ofstream file_out;
    file_out.open(outfilename);
    file_out << quartetOutput;
    file_out.close();

    return (0);
}

