#include "wavelet.h"

void read_signal(double* x,int n,char* filename) {
  int i;
  FILE *fp=fopen(filename,"rt");
  for (i=0;i<n;i++) {
   float t;
   fscanf(fp,"%f\n",&t);
   x[i]=(double)t;
  }
  fclose(fp);
}

void print_signal(double* x,int n) {
  int i;
  for (i=0;i<n;i++) {
    printf("x[%d]=%f\n",i,x[i]);
  }
  printf("\n");
}


void save_signal(double* x,int n,char* filename) {
  int i;
  FILE *fp=fopen(filename,"wt");
  for (i=0;i<n;i++) {
   fprintf(fp,"%f\n",x[i]);
  }
  fclose(fp);
}

double errorIndicator(double* inputSignal, double* outputSignal,int n){
	double error = 0;
	for(int i = 0; i<n; i++){
		error += (inputSignal[i] - outputSignal[i]) * (inputSignal[i] - outputSignal[i]);
	}
	return error;
}

double minArm(double* x, int n){
	double ret = x[0];
	for(int i = 0; i<n; i++){
		if(x[i]<ret){
			ret = x[i];
		}
	}
	return ret;
}


double maxArm(double* x, int n){
	double ret = x[0];
	for(int i = 0; i<n; i++){
		if(x[i]>ret){
			ret = x[i];
		}
	}
	return ret;
}

double meanArm(double* x, int n){
	double ret = 0;
	for(int i = 0; i<n; i++){
		ret += x[i];
	}
	return ret/n;
}


double meanArm2D(double* m ,int xMin, int yMin, int p, int largeur){
	double ret = 0;
	for(int x = xMin; x<p+xMin; x++){
		for(int y = yMin; y<p+yMin; y++){
			ret += m[x+y*largeur];
		}
	}
	return ret/(p*p);
}

double varianceArm2D(double* m ,int xMin, int yMin, int p, int largeur, double mean){
	double ret = 0;
	for(int x = xMin; x<p+xMin; x++){
		for(int y = yMin; y<p+yMin; y++){
			ret += pow(m[x+y*largeur] - mean,2);
		}
	}
	return ret/(p*p);
}
double varMult(double* m, double p, int niveau){
	double ret = 1;
	int i;
	double size;
	double globalSize = pow(p,2);
	double currentSize = p;
	for(i = 1; i<=niveau; i++){
		currentSize = currentSize/2;	
		size = pow(currentSize,2);	

		double vMean = meanArm2D(m, currentSize, 0, currentSize, p);
		double variance = varianceArm2D(m, currentSize, 0, currentSize, p, vMean); 
		ret *= pow(variance, size/globalSize);
		
		vMean = meanArm2D(m,0, currentSize, currentSize,p);
		variance = varianceArm2D(m, 0, currentSize, currentSize, p, vMean);
		ret *= pow(variance, size/globalSize);

		vMean = meanArm2D(m, currentSize, currentSize, currentSize, p);
		variance = varianceArm2D(m, currentSize, currentSize, currentSize, p, vMean); 
		ret *= pow(variance, size/globalSize);
	}

	i--;
	double vMean = meanArm2D(m, 0, 0, currentSize, p);
	double variance = varianceArm2D(m, 0, 0, currentSize, p, vMean);
	ret *= pow(variance, size/globalSize);

	return ret;
}

double debitPixel(double* m, int p, int niveau, double variance, double debitGlobal){
	//std::cout << "Produit : " << varMult(m,p,niveau) << std::endl;
	return debitGlobal + (log2(variance/varMult(m,p,niveau))/2);
}
//////////////////////////////////////////////////////////////////////////////////


void interpolate2(double* x, int n){
	double* y = (double *)malloc(n*sizeof(double));
	memcpy(y,x,n*sizeof(double));

	for(int i = 0; i< n; i++){
		if( (i&1) == 0 ){
			x[i] = y[i/2];
		}else{
			x[i] = 0;
		}
	}

	free(y);
}

void decimate2(double* x, int n){
	double* y = (double *)malloc(n*sizeof(double));
	memcpy(y,x,n*sizeof(double));

	int pair = 0;

	for(int i=0; i<n; i++){
		x[i] = 0;
	}

	for(int i = 0; i< n; i+=2){
		x[i/2] = y[i];
	}


	free(y);
}

void filter(double* x, const double* h, int p, int q){
	double* y = (double *)calloc(p,sizeof(double));

	for(int i = 0; i<p; i++){
		for(int j=0; j<q; j++){
			int k = j-q/2;
			int t = i-k;

			if(t<0){
				t = -t;
			}

			if(t>=p){
				t = (p - 1) - (t-(p-1));
			}

			y[i]+= h[j]*x[t];

		}
	}

	for(int i = 0; i<p; i++){
		x[i] = y[i];
	}

	free(y);
}
////////////////////////////////////////////////////////////////////

void prediction(double* x, int p, double a){
	double* xCpy = (double *)calloc(p,sizeof(double));
	memcpy(xCpy,x,p*sizeof(double));
	for(int i = 0; 2*i+1 < p; i++){ 
		if(p - (2*i+2) <=0){
			x[2*i+1] = a*xCpy[2*i] + xCpy[2*i+1] + a*xCpy[2*i]; 
		}else{

			x[2*i+1] = a*xCpy[2*i] + xCpy[2*i+1] + a*xCpy[2*i+2]; 
		}
	}
	free(xCpy);
}

void misajour(double* x, int p, double a){
	double* xCpy = (double *)calloc(p,sizeof(double));
	memcpy(xCpy,x,p*sizeof(double));

	for(int i = 0; 2*i<p;i++){

		if(2*i-1 <0){
			x[2*i] = a*xCpy[2*i+1] + xCpy[2*i] + a*xCpy[2*i+1];
		}else{
			x[2*i] = a*xCpy[2*i-1] + xCpy[2*i] + a*xCpy[2*i+1];
		}
	}
	free(xCpy);
}

void misalechelle(double* x, int p, double a){
	for(int i = 0; 2*i <p; i++){
		if(p-(2*i+1)>0){
			x[2*i+1] *= a;
		}
		x[2*i] *= 1/a;
	}
}


////////////////////////////////////////////////////////////////////
void analyse(double* x,int p,double* _h0,double* _h1, int filterSize){
	double* xb = (double *)calloc(p,sizeof(double));
	memcpy(xb,x,p*sizeof(double));

	double* xh = (double *)calloc(p,sizeof(double));
	memcpy(xh,x,p*sizeof(double));

	filter(xb,_h0,p,filterSize);
	filter(xh,_h1,p,filterSize);

	decimate2(xb,p);
	decimate2(xh,p);

	for(int i = 0; i<p/2 ; i++){
		x[i] = xb[i];
	}

	for(int i = p/2; i<p ; i++){
		x[i] = xh[i - p/2];
	}

	free(xb);
	free(xh);
}

void analyse_haar(double* x,int p){
	int filterSize = 3;
	double _h0[filterSize] = {1/sqrt(2) , 1/sqrt(2) , 0};
	double _h1[filterSize] = {1/sqrt(2) , -1/sqrt(2) , 0};

	analyse(x,p,_h0,_h1,filterSize);
}

void analyse_97(double* x,int p){
	int filterSize = 9;

	double _h0[filterSize];
	// Filtre biorthogonal 9/7 _h0 (longueur 9)
	_h0[0]=0.037828455507;
	_h0[1]=-0.023849465019;
	_h0[2]=-0.110624404418;
	_h0[3]=0.377402855613;
	_h0[4]=0.852698679009;
	_h0[5]=0.377402855613;
	_h0[6]=-0.110624404418;
	_h0[7]=-0.023849465019;
	_h0[8]=0.037828455507;

	double _h1[filterSize];
	// Filtre biorthogonal 9/7 _h1 (longueur 9)
	_h1[0]=0.064538882629;
	_h1[1]=-0.040689417610;
	_h1[2]=-0.418092273222;
	_h1[3]=0.788485616406;
	_h1[4]=-0.418092273222;
	_h1[5]=-0.040689417610;
	_h1[6]=0.064538882629;
	_h1[7]=0.000000000000;
	_h1[8]=-0.000000000000;

	analyse(x,p,_h0,_h1,filterSize);
}

void analyse_97_lifting(double* x, int p){

	prediction(x,p,-1.586134342);
	misajour(x,p,-0.05298011854);
	prediction(x,p,0.8829110762);
	misajour(x,p,0.4435068522);
	misalechelle(x,p,1/1.149604398);


	double* xCpy = (double *)calloc(p,sizeof(double));
	memcpy(xCpy,x,p*sizeof(double));
	int j = 0;
	for(int i = 0; 2*i<p; i++){
		x[j] = xCpy[2*i];
		j++;
	}
	for(int i = 0; 2*i+1 <p; i++){
		x[j] = xCpy[2*i+1];
		j++;
	}
	free(xCpy);
}

void amr(double* x, int p, int niveau){// 1+j
	if(niveau == 0){
		return;
	}

	analyse_97_lifting(x,p);
	std::cout << "Min de l'approximation de " << niveau << " : " << minArm(x,p/2) <<std::endl;
	std::cout << "Max de l'approximation de " << niveau << " : " << maxArm(x,p/2) <<std::endl;
	std::cout << "Moyenne de l'approximation de " << niveau << " : " << meanArm(x,p/2) <<std::endl;



	std::cout << "Min du détail de " << niveau << " : " << minArm(&(x[p/2]),p/2) <<std::endl;
	std::cout << "Max du détail de " << niveau << " : " << maxArm(&(x[p/2]),p/2) <<std::endl;
	std::cout << "Moyenne du détail de " << niveau << " : " << meanArm(&(x[p/2]),p/2) << "\n" <<std::endl;
	amr(x,p/2,niveau-1);
}

void analyse2D_97(double* m, int p){

	double* mCpy = (double *)calloc(p,sizeof(double));

	for(int y = 0; y<p; y++){
		analyse_97_lifting(&(m[y*p]),p);
	}

	for(int x = 0; x<p; x++){
		for(int y = 0; y< p; y++){
			mCpy[y] = m[x+y*p];
		}
		analyse_97_lifting(mCpy,p);
		for(int y = 0; y<p; y++){
			m[x+y*p] = mCpy[y];
		}
	}

	for(int x = 0; x<p; x++){
		for(int y = 0; y<p; y++){
			//m[x+y*p]*=10;
		}
	}

	free(mCpy);
}

void amr2D_97(double* m,int p,int niveau){// 1+j*3
	double* mCpy = (double *)calloc(p,sizeof(double));
	int i;
	double currentSize = p;
	for(i = 0; i<niveau; i++){
		for(int y = 0; y<p/pow(2,i); y++){
			analyse_97_lifting(&(m[y*p]),p/pow(2,i));
		}

		for(int x = 0; x<p/pow(2,i); x++){
			for(int y = 0; y< p/pow(2,i); y++){
				mCpy[y] = m[x+y*p];
			}
			analyse_97_lifting(mCpy,p/pow(2,i));
			for(int y = 0; y<p/pow(2,i); y++){
				m[x+y*p] = mCpy[y];
			}
		}
		currentSize = currentSize/2;

 		double vMean = meanArm2D(m,currentSize,0,currentSize,p);
 		double variance = varianceArm2D(m,currentSize,0,currentSize,p,vMean);
		//std::cout << "Moyenne de detail horizontal de " << i << " : " << vMean <<std::endl;
		//std::cout << "Variance de detail horizontal de " << i << " : " << variance <<std::endl;
		/*std::cout << "debitPixel de detail horizontal de " << i << " avec b = 2.0: " << debitPixel(m,p,niveau,variance,2.0) <<std::endl;
		std::cout << "debitPixel de detail horizontal de " << i << " avec b = 0.5: " << debitPixel(m,p,niveau,variance,0.5) <<std::endl;
		std::cout << "debitPixel de detail horizontal de " << i << " avec b = 0.25: " << debitPixel(m,p,niveau,variance,0.25) <<std::endl;*/
		std::cout << std::endl;

 		vMean = meanArm2D(m,0,currentSize,currentSize,p);
 		variance = varianceArm2D(m,0,currentSize,currentSize,p,vMean);
		//std::cout << "Moyenne de detail vertical de " << i << " : " << vMean <<std::endl;
		//std::cout << "Variance de detail vertical de " << i << " : " << variance <<std::endl;
		/*std::cout << "debitPixel de detail vertical de " << i << " avec b = 2.0: " << debitPixel(m,p,niveau,variance,2.0) <<std::endl;
		std::cout << "debitPixel de detail vertical de " << i << " avec b = 0.5: " << debitPixel(m,p,niveau,variance,0.5) <<std::endl;
		std::cout << "debitPixel de detail vertical de " << i << " avec b = 0.25: " << debitPixel(m,p,niveau,variance,0.25) <<std::endl;*/
		std::cout << std::endl;


 		vMean = meanArm2D(m,currentSize,currentSize,currentSize,p);
 		variance = varianceArm2D(m,currentSize,currentSize,currentSize,p,vMean);
		//std::cout << "Moyenne de detail diagonal de " << i << " : " << vMean <<std::endl;
		//std::cout << "Variance de detail diagonal de " << i << " : " << variance <<std::endl;
		/*std::cout << "debitPixel de detail diagonal de " << i << " avec b = 2.0: " << debitPixel(m,p,niveau,variance,2.0) <<std::endl;
		std::cout << "debitPixel de detail diagonal de " << i << " avec b = 0.5: " << debitPixel(m,p,niveau,variance,0.5) <<std::endl;
		std::cout << "debitPixel de detail diagonal de " << i << " avec b = 0.25: " << debitPixel(m,p,niveau,variance,0.25) <<std::endl;*/
		std::cout << "------------------------"<< std::endl;
		std::cout << std::endl;

	}

	double vMean = meanArm2D(m,0,0,currentSize,p);
	double variance = varianceArm2D(m,0,0,currentSize,p,vMean);
	//std::cout << "Moyenne de l'approximation de " << i << " : " << vMean <<std::endl;
	//std::cout << "Variance de l'approximation de " << i  << " : " << variance <<std::endl;
	std::cout << "debitPixel de l'approximation de " << i  << " avec b = 2.0: " << debitPixel(m,p,niveau,variance,2.0) <<std::endl;
	std::cout << "debitPixel de l'approximation de " << i  << " avec b = 0.5: " << debitPixel(m,p,niveau,variance,0.5) <<std::endl;
	std::cout << "debitPixel de l'approximation de " << i << " avec b = 0.25: " << debitPixel(m,p,niveau,variance,0.25) <<std::endl;

	free(mCpy);

	for(int x = 0; x<p; x++){
		for(int y = 0; y<p; y++){
			//m[x+y*p]*=5;
		}
	}
}

//////////////////////////////////////////////////////////////////
void synthese(double* x,int p,double* _g0,double* _g1, int filterSize1, int filterSize2){
	double* xbd = (double *)calloc(p,sizeof(double));
	memcpy(xbd,x,(p/2)*sizeof(double));

	double* xhd = (double *)calloc(p,sizeof(double));
	memcpy(xhd,&(x[p/2]),(p/2)*sizeof(double));

	interpolate2(xbd,p);
	interpolate2(xhd,p);

	filter(xbd,_g0,p,filterSize1);
	filter(xhd,_g1,p,filterSize2);

	for(int i = 0; i<p; i++){
		x[i] = xbd[i] + xhd[i];
	}

	free(xbd);
	free(xhd);
}

void synthese_haar(double* x,int p){
	int filterSize = 3;
	double _g0[filterSize] = {0 , 1/sqrt(2), 1/sqrt(2)};
	double _g1[filterSize] = {0 , -1/sqrt(2), 1/sqrt(2)};

	synthese(x,p,_g0,_g1,filterSize,filterSize);
}

void synthese_97(double* x, int p){
	int filterSize1 = 7;
	double _g0[filterSize1];
	// Filtre biorthogonal 9/7 _g0 (longueur 7)
	_g0[0]=-0.064538882629;
	_g0[1]=-0.040689417610;
	_g0[2]=0.418092273222;
	_g0[3]=0.788485616406;
	_g0[4]=0.418092273222;
	_g0[5]=-0.040689417610;
	_g0[6]=-0.064538882629;

	int filterSize2 = 11;
	double _g1[filterSize2];
	// Filtre biorthogonal 9/7 _g1 (longueur 11)
	_g1[0]=0.000000000000;
	_g1[1]=-0.000000000000;
	_g1[2]=0.037828455507;
	_g1[3]=0.023849465019;
	_g1[4]=-0.110624404418;
	_g1[5]=-0.377402855613;
	_g1[6]=0.852698679009;
	_g1[7]=-0.377402855613;
	_g1[8]=-0.110624404418;
	_g1[9]=0.023849465019;
	_g1[10]=0.037828455507;


	synthese(x,p,_g0,_g1,filterSize1,filterSize2);
}

void synthese_97_lifting(double* x, int p){
	double* xCpy = (double *)calloc(p,sizeof(double));
	memcpy(xCpy,x,p*sizeof(double));

	int j = 0;
	for(int i = 0; 2*i<p; i++){
		x[2*i] = xCpy[j];
		j++;
	}
	for(int i = 0; 2*i+1 <p; i++){
		x[2*i+1] = xCpy[j];
		j++;
	}
	free(xCpy);

	misalechelle(x,p,1.149604398);
	misajour(x,p,-0.4435068522);
	prediction(x,p,-0.8829110762);
	misajour(x,p,0.05298011854);
	prediction(x,p,1.586134342);
}


void iamr(double* x, int p, int niveau){

	if(niveau == 0){
		return;
	}

	/*if(niveau == 1){
		synthese_97(x,p);
		return;
	}*/

	synthese_97_lifting(x,p/pow(2,niveau-1));
	iamr(x,p,niveau-1);
}

void iamr2D_97(double* m,int p,int niveau){
	double* mCpy = (double *)calloc(p,sizeof(double));
	for(int i = niveau-1; i>=0; i--){

		for(int x = 0; x<p/pow(2,i); x++){
			for(int y = 0; y< p/pow(2,i); y++){
				mCpy[y] = m[x+y*p];
			}
			synthese_97_lifting(mCpy,p/pow(2,i));
			for(int y = 0; y<p/pow(2,i); y++){
				m[x+y*p] = mCpy[y];
			}
		}

		for(int y = 0; y<p/pow(2,i); y++){
			synthese_97_lifting(&(m[y*p]),p/pow(2,i));
		}
	}
	free(mCpy);
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
int main(){

//////////////////////Decommentez le signal que vous voulez utiliser puis decommenter l'action que vous voulez faire//////////////////////////////	
	/*int n = 4096;
	double* x = (double *)malloc(n*sizeof(double));
	read_signal(x,n,"leleccum.txt");*/

	/*int n = 256;
	double* x = (double *)malloc(n*sizeof(double));
	for(int i = 0; i<n; i++){
		x[i] = i;
	}*/

	/*int n = 512;
	double* x = (double *)malloc(n*sizeof(double));
	read_signal(x,n,"test.txt");*/


//////////////////////Banc de filtres biorthogonaux 9/7//////////////////////////
/*	analyse_97(x,n);
	print_signal(x,n);
	//save_signal(x,n,"analyse97_rampe.txt");
	//save_signal(x,n,"analyse97_lelecum.txt");
	//save_signal(x,n,"analyse97_test.txt");



	synthese_97(x,n);
	print_signal(x,n);
	//save_signal(x,n,"synthese97_rampe.txt");
	//save_signal(x,n,"synthese97_lelecum.txt");
	//save_signal(x,n,"synthese97_test.txt");
*/
/////////////////////////Banc de filtres Haar/////////////////////////////
	/*analyse_haar(x,n);
	print_signal(x,n);
	//save_signal(x,n,"analyseHaar_rampe.txt");
	//save_signal(x,n,"analyseHaar_lelecum.txt");
	//save_signal(x,n,"analyseHaar_test.txt");



	synthese_haar(x,n);
	print_signal(x,n);
	//save_signal(x,n,"syntheseHaar_rampe.txt");
	//save_signal(x,n,"syntheseHaar_lelecum.txt");
	//save_signal(x,n,"syntheseHaar_test.txt");*/

///////////////////////////Lifting biorthogonaux 9/7////////////////////////////////
	/*analyse_97_lifting(x,n);
	print_signal(x,n);
	save_signal(x,n,"analyseLifting_rampe.txt");
	//save_signal(x,n,"analyseLifting_test.txt");
	//save_signal(x,n,"analyseLifting_lelecum.txt");


	synthese_97_lifting(x,n);
	print_signal(x,n);
	save_signal(x,n,"syntheseLifting_lelecum.txt");
	//save_signal(x,n,"syntheseLifting_rampe.txt");
	//save_signal(x,n,"syntheseLifting_test.txt");*/

/////////////////////////////AMR 1D/////////////////////////////////////
/*
	amr(x,n,9);
	//print_signal(x,n);
	save_signal(x,n,"amr9_test.txt");


	iamr(x,n,9);
	//print_signal(x,n);
	save_signal(x,n,"iamr9_test.txt");*/

/////////////////////////////AMR 2D////////////////////////////////////////////////
	/*uint32_t pHauteur,pLargeur;
	double* m = charge_bmp256("./lena.bmp",&pHauteur,&pLargeur);
	int largeur = pLargeur;
	int hauteur = pHauteur;
	amr2D_97(m,hauteur,3);
	//ecrit_bmp256("./lena_analyse_2D_lvl1.bmp",largeur,hauteur,m);
	iamr2D_97(m,hauteur,3);
	//ecrit_bmp256("./lena_Synthese_2D_lvl1.bmp",largeur,hauteur,m);*/

	return 0;
}
