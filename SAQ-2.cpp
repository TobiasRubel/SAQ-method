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
    if (argv[3] == NULL)
    {
        filter = -1;
    }
    else
    {
        filter = stod(argv[3]);
    }

    // READ DATA
    Alignment align = readFASTA(filename);
    map<string, int> tensor = Getcolumns(align);
    int length = align.seq_len;

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
    sp0 = align.taxa[0];
    sp1 = align.taxa[1];
    sp2 = align.taxa[2];
    sp3 = align.taxa[3];

    cout << endl;
    cout << "Taxa 1: " << sp0 << endl;
    cout << "Taxa 2: " << sp1 << endl;
    cout << "Taxa 3: " << sp2 << endl;
    cout << "Taxa 4: " << sp3 << endl;
    cout << "Topologies : \t" << sp0 << sp1 << "|" << sp2 << sp3 << "\t \t "
         << sp0 << sp2 << "|" << sp1 << sp3 << "\t \t "
         << sp0 << sp3 << "|" << sp1 << sp2 << endl;
    cout << "SAQ weights: \t" << score_01 / sum << "\t " << score_02 / sum << "\t " << score_03 / sum << endl;

    ofstream file_out;

    file_out.open(outfilename, std::ios_base::app);

    if (score_01 >= score_02 && score_01 >= score_03)
    {
        cout << "Best: " << sp0 << sp1 << "|" << sp2 << sp3 << endl;
        file_out << sp0 << sp1 << "|" << sp2 << sp3 << " ";
    }
    else if (score_02 >= score_01 && score_02 >= score_03)
    {
        cout << "Best: " << sp0 << sp2 << "|" << sp1 << sp3 << endl;
        file_out << sp0 << sp2 << "|" << sp1 << sp3 << " ";
    }
    else if (score_03 >= score_01 && score_03 >= score_01)
    {
        cout << "Best: " << sp0 << sp3 << "|" << sp1 << sp2 << endl;
        file_out << sp0 << sp3 << "|" << sp1 << sp2 << " ";
    }

    return (0);
}