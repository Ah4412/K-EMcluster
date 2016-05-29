#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace std;

struct vect{
	vector<double> v;

	vect operator+(const vect a) const
	{
		if(a.v.size() == v.size() ){
			vect out;
			for(unsigned i=0; i < a.v.size(); i++){
				out.v.push_back(v[i] + a.v[i]);
			}
			return out;
		}else{
			fprintf(stderr, "Error unmatching vectors '+'\n");
            exit(1);
		}
	}
	vect operator-(const vect a) const
	{
		if(a.v.size() == v.size() ){
			vect out;
			for(unsigned i=0; i < a.v.size(); i++){
				out.v.push_back(v[i] - a.v[i]);
			}
			return out;
		}else{
			fprintf(stderr, "Error unmatching vectors '+'\n");
            exit(1);
		}
	}
};



int main(int argc, char *argv[] ){
	vect nbr1;
	vect nbr2;
	for(unsigned i=0; i<4; i++){
		nbr1.v.push_back(i);
		nbr2.v.push_back(i);
	}
	nbr1 = nbr1 + nbr2;

	for(unsigned i=0; i < nbr1.v.size(); i++)
		std::cout << "Space : " << i << " value : " << nbr1.v[i] << std::endl;

	nbr1 = nbr1 - nbr2;

	for(unsigned i=0; i < nbr1.v.size(); i++)
		std::cout << "Space : " << i << " value : " << nbr1.v[i] << std::endl;

	return 0;
}