//
//  main.cpp
//  hilbert
//
//  Created by Yourui Guo on 2019-09-24.
//  Copyright © 2019 Yourui Guo. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
using namespace std;
#define BITS 5

#define DIM 2
#define ALLBITS BITS*DIM

vector<unsigned int> HilbertAxes(vector<unsigned int> transposedIndex, int dim) {
    vector<unsigned int>X = transposedIndex;
    int n = dim; // n: Number of dimensions
    unsigned int N = 2U << (BITS - 1), P, Q, t;
    int i;
    // Gray decode by H ^ (H/2)
    t = X[n - 1] >> 1;
    // Corrected error in Skilling's paper on the following line. The appendix had i >= 0 leading to negative array index.
    for (i = n - 1; i > 0; i--)
        X[i] ^= X[i - 1];
    X[0] ^= t;
    // Undo excess work
    for (Q = 2; Q != N; Q <<= 1)
    {
        P = Q - 1;
        for (i = n - 1; i >= 0; i--)
            if ((X[i] & Q) != 0U)
                X[0] ^= P; // invert
            else
            {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
    } // exchange
    return X;
}

unsigned int* hilbert(unsigned int hilbertAxes[], int dim){
    unsigned int *X = hilbertAxes;
    int n = dim; // n: Number of dimensions
    
    unsigned int M = 1U << (BITS - 1), P, Q, t;
    int i;
    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1)
    {
        P = Q - 1;
        for (i = 0; i < n; i++)
            if ((X[i] & Q) != 0)
                X[0] ^= P; // invert
            else
            {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
    } // exchange
    // Gray encode
    for (i = 1; i < n; i++)
        X[i] ^= X[i - 1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if ((X[n - 1] & Q)!=0)
            t ^= Q - 1;
    for (i = 0; i < n; i++)
        X[i] ^= t;
    
    return X;
}



int hilbertIndexToAxes(unsigned int hilbertAxes[], int dim){
    stringstream ss;
    string axes[dim];
    std::string bi;
    for (int i=0; i<dim; i++) {
        stringstream().swap(ss);
        ss << std::bitset<ALLBITS>(hilbertAxes[i]);
        axes[i] = ss.str();
    }
    //cout<<axes[0] << ' ' <<axes[1]<< endl;
    for (int i=0; i<BITS; i++) {
        for (int j=0; j<dim; j++) {
            bi+=axes[j][i];
        }
    }
    return int(std::stoull(bi, 0, 2));
    
}

vector<unsigned int> hilbertAxesToIndex(int index, int dim) {
    char axes[dim][BITS+1];
    vector<unsigned int>ind;
    string s =  std::bitset<ALLBITS>(index).to_string();
    int a = int(s.size()-1);

    for (int i=BITS-1; i>=0; i--) {
        for (int j=dim-1; j>=0; j--) {
            axes[j][i] = s[a--];
        }
    }
    
    for (int i=0; i<dim; i++) {
        axes[i][BITS] = '\0';
        ind.push_back(uint(std::stoull(axes[i], 0, 2))) ;
        //cout<<std::stoull(axes[i], 0, 2)<<endl;
        
    }
    return ind;
}

int main(int argc, const char * argv[]) {
    unsigned int axes[2];
    vector<unsigned int> temp, h;

    //axes[0] = 2;
    //axes[1] = 1;
    //unsigned int *t = hilbert(axes, 2);
    //int index = hilbertIndexToAxes(t, 2);
    for (int i=0; i<20; i++) {
        int index = i;
        temp = hilbertAxesToIndex(index, 2);
        h = HilbertAxes(temp, 2);
        cout <<h[0] <<' ' <<h[1]<<endl;
    }
    
    
    return 0;
}
