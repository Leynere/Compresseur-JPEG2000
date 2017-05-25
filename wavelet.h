#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include "bmp.c"
#include "quantlm.cpp"

void save_signal(double* x,int n,char* filename);
void print_signal(double* x,int n);		
void read_signal(double* x,int n,char* filename);
double errorIndicator(double* inputSignal, double* outputSignal, int n);

void interpolate2(double* x, int n);
void decimate2(double* x, int n);
void filter(double* x, const double* h, int p, int q);

void analyse(double* x,int p,double* _h0,double* _h1, int filterSize);
void analyse_haar(double* x,int p);
void analyse_97(double* x,int p);

void synthese(double* x,int p,double* _g0,double* _g1, int filterSize1, int filterSize2);
void synthese_haar(double* x,int p);
void synthese_97(double* x, int p);

void analyse_97_lifting(double* x,int p);
void synthese_97_lifting(double* x, int p);

void amr(double* x, int p, int niveau);
void iamr(double* x, int p, int niveau);

void analyse2D_97(double* signal, int p);
void amr2D_97(double* m,int p,int j);