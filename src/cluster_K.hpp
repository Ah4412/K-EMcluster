#ifndef cluster_K_hpp
#define cluster_K_hpp

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

void cluster_K(const int nbr_cls, const int size, const int range, coordonnees_p *p){
	coordonnees_c c[nbr_cls];


	int quotient = (int)size/nbr_cls;
	int rest = size%nbr_cls;
	int r2 = (int)((range*10)/100);

	for(unsigned i=0; i<nbr_cls; i++){
		c[i].x = rand() % range;
		c[i].y = rand() % range;

		for(unsigned v=i*quotient; v<(i+1)*quotient; v++){
			p[v].x= c[i].x + rand()%r2 - (int)r2/2;
			p[v].y= c[i].y + rand()%r2 - (int)r2/2;
		}
	} // Random clustered created.

	for(unsigned i= size-rest; i<size; i++){
		p[i].x = rand() % range;
		p[i].y = rand() % range;
	} // set random points free from clustered as noise 
}	  // if size%nbr_cls != 0

#endif