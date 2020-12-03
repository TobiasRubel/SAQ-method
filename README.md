# SAQ: semi-algebraic quartet reconstruction method

SAQ V.1.0
Update: Dec/2020

This code has been written by M. Garrote-López. 
Please quote the paper: M. Casanellas, J. Fernández-Sánchez, M. Garrote-López, "SAQ: semi-algebraic quartet reconstruction method", https://arxiv.org/abs/2011.13968

## Description
The SAQ method has been implemented in C++ to be run in a UNIX terminal. Given an alignment of four sequences, SAQ computes three weights, one for each possible tree topology. 


### Input

The input of the method is a FASTA file containing the alignment of four species and the parameter _filter_ which is optional.


| Input | Description | Example |
|:------------- |:------------- | :-----: |
| fastaFile | Path to the fasta file with the alignment | /home/user/data/alignment.fa |
| filter (optional)| Minimum value accepted on the entries of the vectors obtained by a _T_-leaf transformation on the distribution obtained from the alignment. <br /> The default value is _filter = -1_ | -0.5 |


### Output

The program outputs three weights provided by the method SAQ applied to the alignment. Each weight is computed according to one possible topology relating the sequences of the alignment _seq1, seq2, seq3, seq4_:

* Weight of tree _seq1,seq2 | seq3,seq4_
* Weight of tree _seq1,seq3 | seq2,seq4_
* Weight of tree _seq1,seq4 | seq2,seq3_ 


## Usage 

The executable files SAQ have been written in C++ code and compiled using the c++11 compiler. To run SAQ it is necessary to have installed the C++ linear algebra library Armadillo (http://arma.sourceforge.net/). 

### Compiling SAQ
```
  g++ -std=c++11 -c functionsSAQ.cpp -O1 -larmadillo
  g++ -std=c++11 SAQ.cpp functionsSAQ.o -o SAQ -O1 -larmadillo
```
    
### Running an experiment
```
  ./SAQ fastaFile filter
```
The parameter _filter_ can be omitted. 

### Examples:
```
  ./SAQ /home/user/data/alignment.fa -0.5
  ./SAQ /home/user/data/alignment.fa
```
