#ifndef K_hpp
#define K_hpp

template<class T>
void genrdm(T& p, int N){
	for(unsigned i=0; i< N; i++){
		p[i].x = rand() % 1000;
		p[i].y = rand() % 1000;
	}
	return;
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
	return; // if size%nbr_cls != 0
}	

template<class Tp, class Tc>
void closest_c(Tp& p, Tc& c, int numb_cluster){
	double dx, dy, dist;
	double dmin;
	int index;

	for(unsigned v=0; v<numb_cluster; v++){
		c[v].n = 0;
	}

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
		c[index].n++;
	}
	return;
}

template<class Tp, class Tc>
void recalculate_c(Tp& p, Tc& c, int numb_cluster){
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
	return;
}

template<class Tp, class Tc>
double return_J(Tp& p, Tc& c, int numb_cluster){
	double J=0;
	double dx=0, dy=0;
	for(unsigned i=0; i<table_size; i++){
		for(unsigned v=0; v<numb_cluster; v++){
			if(p[i].rc[v] == 1){
				dx = (p[i].x - c[v].x)*(p[i].x - c[v].x);
				dy = (p[i].y - c[v].y)*(p[i].y - c[v].y);
				J += dx + dy;
			}
		}
	}
	return J;
}

double fixclus(coordonnees_p *p, const int nbr_cls, vector<coordonnees_c> &C){
	double J, Jtmp, Jfinal;
	bool first = true;
	
	vector<coordonnees_c> cluster(nbr_cls); 
	vector<coordonnees_c> tmp(nbr_cls);
	
	for(unsigned t=0; t<10; t++){

		/////Initialize starting set (do it 10 times)
		bool J_smaller = true;
		genrdm(cluster, nbr_cls);
		closest_c(p, cluster, nbr_cls);
		J = return_J(p, cluster, nbr_cls);

			while(J_smaller){
				J_smaller = false;

				closest_c(p, cluster, nbr_cls);
				recalculate_c(p, tmp, nbr_cls);
				Jtmp = return_J(p, tmp, nbr_cls);

				if(Jtmp < J){
					J = Jtmp;
					J_smaller = true;
					swap(cluster, tmp); 
				}
			}

			if((first == true) || (J < Jfinal)){
				first = false;
				Jfinal = J;
				swap(C, cluster);
			}
	}

	return Jfinal;
}

double return_AIC(int k, double J){ //{AIC} =2q(k)-2\ln(L) q(k) = Dk
	double AIC = 0.4*(double)k + log(J);
	return AIC;
}

void AIC_K(int &nbr_cls_st, vector<coordonnees_c> &C, coordonnees_p *p){

	double J, Jtmp;
	double AIC, AICtmp;

	coordonnees_c x;
	x.x = 0;
	x.y = 0;

	vector<coordonnees_c> Ctmp(nbr_cls_st);

	/////Do K mean for initial number of cluster.
	J = fixclus(p, nbr_cls_st, C);
	AIC = return_AIC(nbr_cls_st, J);
	//cout << "Value of J : "<< J << " and log : "<< log(J) << endl;
	//cout << "Value of AIC : "<< AIC << endl << endl;
	/////Loop incrementing number of cluster and checking the AIC.
	bool AIC_smaller = true;

	while(AIC_smaller){
		AIC_smaller = false;

		nbr_cls_st++;
		Ctmp.push_back(x);

		for(unsigned i = 0; i < table_size; i++){
			p[i].rc.push_back(0);
		}

		Jtmp = fixclus(p, nbr_cls_st, Ctmp);
		AICtmp = return_AIC(nbr_cls_st, Jtmp);	
		//cout << "Value of J : "<< Jtmp << " and log : "<< log(J) << endl;
		//cout << "Value of AIC : "<< AICtmp << endl << endl;


		if(AICtmp < AIC){
			J = Jtmp;
			AIC = AICtmp;
			AIC_smaller = true;
			C.push_back(x);
			swap(C, Ctmp);
		}else{
			nbr_cls_st--;
			for(unsigned i=0; i<table_size; i++)
				p[i].rc.pop_back();
			closest_c(p, C, nbr_cls_st);
			return;
		}
	}
}

#endif