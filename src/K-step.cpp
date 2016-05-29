#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <python2.7/Python.h>

#define table_size 1000
#define numb_cluster 6

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




int main(int argc, char *argv[]){
	cout << "program started" << endl; //flush

	srand (time(NULL));
	double J, Jtmp;
	bool J_smaller = true;
	PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue, *pXVec, *pYVec;
    PyObject *pX, *pY;
    int i=0;

    const char *scriptDirectoryName = "/home/rishie/Documents/Project/src/python/";


	coordonnees_p *p = new coordonnees_p[table_size];
	if(p == NULL){
		cout << "couldn't allocate p" << endl;
		return 0;
	}
	coordonnees_c cluster[numb_cluster]; // stack problem related (limited)
	coordonnees_c tmp[numb_cluster];

	cluster_K(numb_cluster, table_size, 1000, p);
	genrdm(cluster, numb_cluster);

	cout << "created random set" << endl;

	closest_c(p, cluster);
	J = return_J(p, cluster);


	Py_Initialize();

    PyObject *sysPath = PySys_GetObject((char *)"path");
    PyObject *path = PyString_FromString(scriptDirectoryName);
    int result = PyList_Insert(sysPath, 0, path);
    pModule = PyImport_ImportModule("plot");

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, "plot");
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {

            pXVec = PyTuple_New(numb_cluster); 
            pYVec = PyTuple_New(numb_cluster);
            pX = PyTuple_New(table_size);
            pY = PyTuple_New(table_size);

            for (i = 0; i < table_size; i++) {
                pValue = PyFloat_FromDouble(p[i].x);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument x\n");
                    return 0;
                }
                PyTuple_SetItem(pX, i, pValue);
                //set Y's value on other array.
                pValue = PyFloat_FromDouble(p[i].y);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument y\n");
                    return 0;
                }
                PyTuple_SetItem(pY, i, pValue);   
            }

            	while(J_smaller){
					J_smaller = false;

					 for (i = 0; i < numb_cluster; i++) {
			                pValue = PyFloat_FromDouble(cluster[i].x);
			                if (!pValue) {
			                    Py_DECREF(pArgs);
			                    Py_DECREF(pModule);
			                    fprintf(stderr, "Cannot convert argument x\n");
			                    return 0;
			                }
			                PyTuple_SetItem(pXVec, i, pValue);
			                //set Y's value on other array.
			                pValue = PyFloat_FromDouble(cluster[i].y);
			                if (!pValue) {
			                    Py_DECREF(pArgs);
			                    Py_DECREF(pModule);
			                    fprintf(stderr, "Cannot convert argument y\n");
			                    return 0;
			                }
			                PyTuple_SetItem(pYVec, i, pValue);   
			            }

			            pArgs = PyTuple_New(4);

			            PyTuple_SetItem(pArgs, 0, pXVec);
						PyTuple_SetItem(pArgs, 1, pYVec);
						PyTuple_SetItem(pArgs, 2, pX);
            			PyTuple_SetItem(pArgs, 3, pY);

						pValue = PyObject_CallObject(pFunc, pArgs);

						Py_INCREF(pArgs);

					closest_c(p, cluster);
					recalculate_c(p, tmp);
					Jtmp = return_J(p, tmp);

					if(Jtmp < J){
						J = Jtmp;
						J_smaller = true;
						swap(cluster, tmp); 

						cout <<endl;
						for(unsigned i=0; i<numb_cluster; i++){
							cout << "Cluster nbr : " << i << " -> x = " << cluster[i].x << " y = " << cluster[i].y << endl;
						}
						cout << "Value of J : " << J << endl;
						cout <<endl;
					}	
				}
				Py_DECREF(pArgs);

            if (PyInt_AsLong(pValue) == 1) {
                printf("Result of call, printed.\n");
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 0;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function/file\n");
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load module\n");
        return 0;
    }
    Py_Finalize();



	delete [] p;

	
	/*Py_SetProgramName(argv[0]);
  	Py_Initialize();
 	PyRun_SimpleString("from time import time,ctime\n"
                     "print 'Today is',ctime(time())\n");
 	Py_Finalize();*/

	std::cout << "The end!" << std::endl;
	return 0;
}
	
