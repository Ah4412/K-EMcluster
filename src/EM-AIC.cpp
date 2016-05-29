#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <python2.7/Python.h>

#define PI 2*3.141592653589793//23846
#define table_size 1000

using namespace std;

struct coordonnees_p{
	double x;
	double y;
	vector<double> rc;
};

struct coordonnees_c{
	double x;
	double y;
	double cov[2][2];
	double inv[2][2];
	double det;
	double pi;
	int n;
};

#include "EM_cluster.hpp"
#include "plot_dyn.hpp"
#include "K.hpp"
#include "EM.hpp"
	
bool filter_K_cluster(int &nbr_cls_st, vector<coordonnees_c> &C, coordonnees_p *p, int n=0){
	bool anychange = false;
	for(unsigned v = 0; v< nbr_cls_st; v++){
		if(C[v].n <= 4){		//Percentage of dataa point and expected cluster ??
			C[v] = C.back();
			C.pop_back();
			nbr_cls_st--;
			for(unsigned i=0; i<table_size; i++){
				p[i].rc.pop_back();
			}
			closest_c(p, C, nbr_cls_st);
			anychange = true;
		}
	} 
	return anychange;
}


int main(int argc, char *argv[]){
	if(argc != 3){
		cout << "Unexpected number of arguments : " << argc << " stopping program." << endl;
		exit(0);	
	}
	srand(time(NULL));
	/////Define variables.
	int nbr_cls_st = atoi(argv[1]);
	int nbr_cls = atoi(argv[2]);

	vector<coordonnees_c> C(nbr_cls_st);
	coordonnees_c x;

	/////Define data set (randomized)
	coordonnees_p *p = new coordonnees_p[table_size];
	if(p == NULL){
		cout << "couldn't allocate p" << endl;
		return 0;
	}
	else{
		cluster_K(nbr_cls, table_size, 1000, p);
		EM_cluster_p(nbr_cls_st, (int)20);
	}

	for(unsigned i = 0; i < table_size; i++){
		for(unsigned v=0; v < nbr_cls_st; v++){
			p[i].rc.push_back(0);
		}
	}

	bool anychange;
	double L, Ltmp;
	double BIC, BICtmp;

	AIC_K(nbr_cls_st, C, p); // initialize around the right number of cluster use fixclus after.
	/*for(unsigned v=0; v < nbr_cls_st; v++){
		cout << "Number of point in cluster " << v <<" : " << C[v].n << endl;
		cout << "co-ordinate x: " << C[v].x << " & y: " << C[v].y <<endl << endl;
	}*/
	anychange = filter_K_cluster(nbr_cls_st, C, p);
	cout << endl << endl;

	L = fix_EM(nbr_cls_st, C, p, true);
	BIC = return_BIC(nbr_cls_st, L);

	cout << "step 1, done." << endl;

	vector<coordonnees_c> Ctmp(nbr_cls_st);
	bool BIC_smaller = true;
	bool upward =false;

	while(BIC_smaller){
		BIC_smaller = false;

		nbr_cls_st++;
		Ctmp.push_back(x);

		for(unsigned i = 0; i < table_size; i++){
			p[i].rc.push_back(0);
		}

		double J = fixclus(p, nbr_cls_st, Ctmp);

		anychange = filter_K_cluster(nbr_cls_st, Ctmp, p);
		if(anychange){
			resp(nbr_cls_st, C, p);
			break;
		}

		Ltmp = fix_EM(nbr_cls_st, Ctmp, p, true);
		BICtmp = return_BIC(nbr_cls_st, Ltmp);	


		if(BICtmp < BIC){
			L = Ltmp;
			BIC = BICtmp;
			BIC_smaller = true;
			upward = true;
			C.push_back(x);
			swap(C, Ctmp);
			cout << "step 2, increment." << endl;
		}else{	
			nbr_cls_st--;
			for(unsigned i=0; i<table_size; i++)
					p[i].rc.pop_back();
			resp(nbr_cls_st, C, p);
			cout << "step 2, done." << endl;
			break;
		}
	}

	if(!upward){
		vector<coordonnees_c> Ctp(nbr_cls_st);
		BIC_smaller = true;

		while(BIC_smaller){
			BIC_smaller = false;

			nbr_cls_st--;
			Ctp.pop_back();

			for(unsigned i = 0; i < table_size; i++){
				p[i].rc.pop_back();
			}

			double J = fixclus(p, nbr_cls_st, Ctp);

			Ltmp = fix_EM(nbr_cls_st, Ctp, p, true);
			BICtmp = return_BIC(nbr_cls_st, Ltmp);	


			if(BICtmp < BIC){
				L = Ltmp;
				BIC = BICtmp;
				BIC_smaller = true;
				upward = true;
				C.pop_back();
				swap(C, Ctp);
				cout << "step 3, decrement." << endl;
			}else{	
				nbr_cls_st++;
				for(unsigned i=0; i<table_size; i++)
						p[i].rc.push_back(0);
				resp(nbr_cls_st, C, p);
				cout << "step 3, done." << endl;
				break;
			}
		}
	}
	cout << "Number of cluster : " << C.size() << endl;
	/////Release the data set as program end.
	plot_dyn_p(C, p, nbr_cls_st);

	delete [] p;

	cout << endl;
	return 0;
}