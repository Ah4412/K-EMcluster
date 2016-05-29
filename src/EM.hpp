#ifndef EM_hpp
#define EM_hpp

#define PI 2*3.141592653589793//23846

bool init_covar(const int nbr_cls_st, vector<coordonnees_c> &C, const coordonnees_p *p){
	// = 1/n * sum( (xn - uk)* (xn -uk)T )over N after K mean assume 
	// resp to be 1 and 0 for others..n is already set.
	int n = 0;

	for(unsigned v=0; v<nbr_cls_st; v++){
		C[v].cov[0][0] = 0;
		C[v].cov[0][1] = 0;
		C[v].cov[1][0] = 0;
		C[v].cov[1][1] = 0;
	}

	for(unsigned i= 0; i<table_size; i++){
		for(unsigned v=0; v<nbr_cls_st; v++){
			if(p[i].rc[v] == 1){
				C[v].cov[0][0] += (p[i].x - C[v].x)*(p[i].x - C[v].x);
				C[v].cov[0][1] += (p[i].x - C[v].x)*(p[i].y - C[v].y);
				C[v].cov[1][0] += (p[i].y - C[v].y)*(p[i].x - C[v].x);
				C[v].cov[1][1] += (p[i].y - C[v].y)*(p[i].y - C[v].y);
				n++;
			}
		}
	}

	//cout << "number of points assessed :" << n <<endl;

	for(unsigned v=0; v<nbr_cls_st; v++){
		if(C[v].n <= 1){
			//cout << "Cluster with no assignement" << endl;
			return false;
		}

		C[v].cov[0][0] /= (double)C[v].n;
		C[v].cov[0][1] /= (double)C[v].n;
		C[v].cov[1][0] /= (double)C[v].n;
		C[v].cov[1][1] /= (double)C[v].n;
		C[v].pi = (double)C[v].n / table_size;
		C[v].det = C[v].cov[0][0]*C[v].cov[1][1] - C[v].cov[1][0]*C[v].cov[0][1];

		if(C[v].det == 0){
			//cout << "Number : " << C[v].n << endl;
			//cout << "Initialization probleme : det = " << C[v].det << endl;
			//cout << "Then : "<< C[v].cov[0][0]<<" "<<C[v].cov[1][1]<<" "<<C[v].cov[1][0]<< " " <<C[v].cov[0][1]<<endl;
			return false;
		}

		C[v].inv[0][0] = C[v].cov[1][1] / C[v].det;
		C[v].inv[1][0] = -1*(C[v].cov[0][1] / C[v].det);
		C[v].inv[0][1] = C[v].inv[1][0];
		C[v].inv[1][1] = C[v].cov[0][0] / C[v].det;
		C[v].det = sqrt(C[v].det);
	}
	return true;
}

double return_G(const coordonnees_p *p, const coordonnees_c &C){
	//3 values to calculate, first the matrix reduction, then exp then the determinant 
	double det = 1/(PI*C.det);
	double a = p->x - C.x;
	double b = p->y - C.y;
	double e = -0.5 * (a*a*C.inv[0][0] + 2*a*b*C.inv[0][1] + b*b*C.inv[1][1]);

	if(det == 0){
		cout << "return_G probleme : sum = 0" << endl;
		exit(1);
	}

	a = det*exp(e);
	return a; //samata (aka most awesome person around was here)
}

double log_l(const int nbr_cls_st, vector<coordonnees_c> &C,  coordonnees_p *p){
	double L = 0;
	for(unsigned i=0; i<table_size; i++){
		double sum =0;
		for(unsigned v=0; v<nbr_cls_st; v++){
			sum += C[v].pi*return_G(&p[i], C[v]);
		}
		if(sum == 0){
			cout << "log_l probleme : sum = 0" << endl;
			exit(1);
		}
		//cout << "Sum : " << sum << " and log : " << log(sum)<<endl;
		L += log(sum);
	}
	return L;
}

void resp(const int nbr_cls_st, vector<coordonnees_c> &C,  coordonnees_p *p){
	// for each point go through each cluster and retrieve gaussian value...
	// Gaussian value : 1/(2pi*root(det) ) * w

	double buff[nbr_cls_st];
	for(unsigned v=0; v<nbr_cls_st;v++){
		buff[v] = 0;
	}

	for(unsigned i =0; i<table_size; i++){
		double sum = 0;
		double G[nbr_cls_st];
		for(unsigned v =0; v < nbr_cls_st; v++){
			G[v] = return_G(&p[i], C[v]);
			sum += C[v].pi*G[v];
			//cout << G[v] << " " << C[v].pi << " "<< C[v].n << " ";
			//cout << p[i].rc[v] << " ";
		}
		//cout << "sum = " << sum << endl;

		for(unsigned v=0; v<nbr_cls_st; v++){
			p[i].rc[v] = (C[v].pi*G[v])/sum;
			buff[v] += p[i].rc[v];
			//cout << p[i].rc[v] << " ";
		}
	}

	for(unsigned v=0; v<nbr_cls_st; v++){
		C[v].pi = buff[v]/ (double)table_size;
		//cout << " " << C[v].pi <<endl;
	}
}

void mean_EM(const int nbr_cls_st, vector<coordonnees_c> &C, const coordonnees_p *p){
	for(unsigned v=0; v<nbr_cls_st; v++){
		C[v].x = 0;
		C[v].y = 0;		
	}

	for(unsigned i=0; i<table_size; i++){
		for(unsigned v=0; v<nbr_cls_st; v++){
			C[v].x += p[i].rc[v]*p[i].x;
			C[v].y += p[i].rc[v]*p[i].y;
		}
	}
	
	for(unsigned v=0; v<nbr_cls_st; v++){
			C[v].x /= (C[v].pi*table_size);
			C[v].y /= (C[v].pi*table_size);
	}
}

void covar_EM(const int nbr_cls_st, vector<coordonnees_c> &C, const coordonnees_p *p){
	// = 1/n * sum( (xn - uk)* (xn -uk)T )over N after K mean assume 
	// resp to be 1 and 0 for others..n is already set.
	for(unsigned v=0; v<nbr_cls_st; v++){
		C[v].cov[0][0] = 0;
		C[v].cov[0][1] = 0;
		C[v].cov[1][0] = 0;
		C[v].cov[1][1] = 0;
	}

	for(unsigned i= 0; i<table_size; i++){
		for(unsigned v=0; v<nbr_cls_st; v++){
			C[v].cov[0][0] += p[i].rc[v]*(p[i].x - C[v].x)*(p[i].x - C[v].x);
			C[v].cov[0][1] += p[i].rc[v]*(p[i].x - C[v].x)*(p[i].y - C[v].y);
			C[v].cov[1][0] += p[i].rc[v]*(p[i].y - C[v].y)*(p[i].x - C[v].x);
			C[v].cov[1][1] += p[i].rc[v]*(p[i].y - C[v].y)*(p[i].y - C[v].y);				
		}
	}

	for(unsigned v=0; v<nbr_cls_st; v++){
		C[v].cov[0][0] /= (C[v].pi*table_size);  //pi*N = Nk 
		C[v].cov[0][1] /= (C[v].pi*table_size);
		C[v].cov[1][0] /= (C[v].pi*table_size);
		C[v].cov[1][1] /= (C[v].pi*table_size);
		//C[v].pi = (double)C[v].n / table_size;
		C[v].det = C[v].cov[0][0]*C[v].cov[1][1] - C[v].cov[1][0]*C[v].cov[0][1];
		C[v].inv[0][0] = C[v].cov[1][1] / C[v].det;
		C[v].inv[1][0] = -1*(C[v].cov[0][1] / C[v].det);
		C[v].inv[0][1] = C[v].inv[1][0];
		C[v].inv[1][1] = C[v].cov[0][0] / C[v].det;
		C[v].det = sqrt(C[v].det);
	}
	return;
}

double fix_EM(const int nbr_cls_st, vector<coordonnees_c> &Cfinal, coordonnees_p *p, bool postK = false){
	double L, Ltmp, Lfinal;
	bool first = true;

	vector<coordonnees_c> C(nbr_cls_st);
	vector<coordonnees_c> Ctmp(nbr_cls_st);

	for(unsigned t=0; t<10; t++){
		bool start =false;

		if(!postK){
			while(!start){
				genrdm(C, nbr_cls_st);
				closest_c(p, C, nbr_cls_st);
				start = init_covar(nbr_cls_st, C, p);
			}
		}else{
			start = init_covar(nbr_cls_st, Cfinal, p);
			C = Cfinal;
			if(!start){
				//cout << "AIC_K start initialization didn't work." << endl << "Jump to normal EM" << endl;
				postK =false;
				while(!start){
					genrdm(C, nbr_cls_st);
					closest_c(p, C, nbr_cls_st);
					start = init_covar(nbr_cls_st, C, p);
				}
			}
		}

		L = log_l(nbr_cls_st, C, p);
		Lfinal = L;
		
		bool L_smaller = true;
		Ctmp = C;
		int n=0;

		while(L_smaller){
			L_smaller = false;

			resp(nbr_cls_st, Ctmp, p);
			mean_EM(nbr_cls_st, Ctmp, p);
			covar_EM(nbr_cls_st, Ctmp, p);
			Ltmp = log_l(nbr_cls_st, Ctmp, p);

			if(Ltmp > L){
				n++;
				L_smaller = true;
				L = Ltmp;
				C = Ctmp;
			}
		}
		if((first == true && !postK) || (L > Lfinal)){
			first = false;
			Lfinal = L;
			Cfinal = C;
		}
		if(postK){
			//cout << "Iteration : " << n << endl;
			break;
		}

	}
	return Lfinal;
}

double return_BIC(const int k, const double L){
	double BIC = -1*L + 4*log(table_size)*k;
	return BIC;
}

void BIC_EM(int &nbr_cls_st, vector<coordonnees_c> &C, coordonnees_p *p){
	double L , Ltmp;
	double BIC, BICtmp;

	coordonnees_c x;
	vector<coordonnees_c> Ctmp(nbr_cls_st);
	
	L = fix_EM(nbr_cls_st, C, p);
	BIC = return_BIC(nbr_cls_st, L);
	/*cout << "Cost : " << L << endl;
	cout << "BIC : " << BIC << endl<< endl;*/


	bool BIC_smaller = true;

	while(BIC_smaller){
		BIC_smaller = false;

		nbr_cls_st++;
		Ctmp.push_back(x);

		for(unsigned i = 0; i < table_size; i++){
			p[i].rc.push_back(0);
		}

		Ltmp = fix_EM(nbr_cls_st, Ctmp, p);
		BICtmp = return_BIC(nbr_cls_st, Ltmp);	
		cout << "Value of L : "<< Ltmp << endl;
		cout << "Cost of CLuster : " << log(table_size)*nbr_cls_st << endl;
		cout << "Value of AIC : "<< BICtmp << endl << endl;


		if(BICtmp < BIC){
			L = Ltmp;
			BIC = BICtmp;
			//if(nbr_cls_st < 7){
				BIC_smaller = true;
			//}
			C.push_back(x);
			swap(C, Ctmp);
		}else{
			nbr_cls_st--;
			for(unsigned i=0; i<table_size; i++)
				p[i].rc.pop_back();
			resp(nbr_cls_st, C, p);
			return;
		}
	}
	return;
}

#endif