#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <python2.7/Python.h>

#define table_size 1000
#define numb_cluster 7

using namespace std;

struct coordonnees_p{
	double x;
	double y;
	int rc[numb_cluster];
};

struct coordonnees_c{
	double x;
	double y;
};

#include "cluster_K.hpp"
#include "plotpy.hpp"
#include "plot.hpp"
/*//////////////////////////////
		template<class T>
		void genrdm(T& p, int N){
			for(unsigned i=0; i< N; i++){
				p[i].x = rand() % 1000;
				p[i].y = rand() % 1000;
			}
		}

		template<class Tp, class Tc>
		void closest_c(Tp& p, Tc& c){
			double dx, dy, dist;
			double dmin;
			int index = 3;

			for(unsigned i=0; i<table_size; i++){
				dmin=10000000;
				for(unsigned v=0; v<numb_cluster; v++){
					dx = (p[i].x - c[v].x)*(p[i].x - c[v].x);
					dy = (p[i].y - c[v].y)*(p[i].y - c[v].y);
					dist = dx + dy;
					if(dist < dmin){
						index = v;
						dmin = dist;
					}
					p[i].rc[v] = 0;
				}
				p[i].rc[index] = 1;
			}
		}

		template<class Tp, class Tc>
		void recalculate_c(Tp& p, Tc& c){
			double avx, avy;
			double count;
			for(unsigned v=0; v<numb_cluster; v++){
				avx =0;
				avy =0;
				count =0;
				for(unsigned i=0; i<table_size; i++){
					if(p[i].rc[v] == 1){
						avx += p[i].x;
						avy += p[i].y;
						count++;
					}
				}
				if(count != 0){
					c[v].x = avx/count;
					c[v].y = avy/count;
				}
			}
		}

		template<class Tp, class Tc>
		double return_J(Tp& p, Tc& c){
			double J=0;
			double dx=0, dy=0;
			for(unsigned i=0; i<table_size; i++){
				for(unsigned v=0; v<numb_cluster;v++){
					if(p[i].rc[v] == 1){
						dx = (p[i].x - c[v].x)*(p[i].x - c[v].x);
						dy = (p[i].y - c[v].y)*(p[i].y - c[v].y);
						J += dx + dy;
					}
				}
			}
			return J;
		}
*/
int main(int argc, char *argv[]){
	cout << "program started" << endl; //flush

	srand (time(NULL));
	double J, Jtmp, Jfinal;
	bool first = true;

	coordonnees_p *p = new coordonnees_p[table_size];
	if(p == NULL){
		cout << "couldn't allocate p" << endl;
		return 0;
	}
	coordonnees_c final[numb_cluster] = {0};
	coordonnees_c cluster[numb_cluster] = {0}; // stack problem related (limited)
	coordonnees_c tmp[numb_cluster] = {0};

	cout << "space allocated" << endl;

	cluster_K(numb_cluster, table_size, 1000, p);

	for(unsigned t=0; t<10; t++){

		bool J_smaller = true;
		genrdm(cluster, numb_cluster);

		cout << "created random set" << t << endl;

		closest_c(p, cluster);
		J = return_J(p, cluster);

			while(J_smaller){
				J_smaller = false;

				closest_c(p, cluster);
				recalculate_c(p, tmp);
				Jtmp = return_J(p, tmp);

				if(Jtmp < J){
					J = Jtmp;
					J_smaller = true;
					swap(cluster, tmp); 
				}	
			}

			for(unsigned i=0; i<numb_cluster; i++){
						cout <<"Cluster nbr : " << i << " -> x = " << cluster[i].x << " y = " << cluster[i].y << endl;
			}
			cout << "Value of J : " << J << endl;
			cout <<endl;

			if((first == true) || (J < Jfinal)){
				first = false;
				Jfinal = J;
				swap(final, cluster);
			}
	}

	cout <<"Final cluster : " << endl;
	for(unsigned i=0; i<numb_cluster; i++){
		cout <<"Cluster nbr : " << i << " -> x = " << final[i].x << " y = " << final[i].y << endl;
	}
	cout << "Value of J :"<<Jfinal<< endl;

	plot_p(final, p);


	delete [] p;

	
	/*Py_SetProgramName(argv[0]);
  	Py_Initialize();
 	PyRun_SimpleString("from time import time,ctime\n"
                     "print 'Today is',ctime(time())\n");
 	Py_Finalize();*/

	std::cout << "The end!" << std::endl;
	return 0;
}
	
