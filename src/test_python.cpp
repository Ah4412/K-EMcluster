#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <python2.7/Python.h>

#define numb_cluster 5

using namespace std;

struct coordonnees_c{
	double x;
	double y;
};

template<class T>
void genrdm(T& p, int N){
	for(unsigned i=0; i< N; i++){
		p[i].x = rand() % 1000;
		p[i].y = rand() % 1000;
	}
}

int main(int argc, char *argv[]){

	PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    int i;


	srand (time(NULL));

	coordonnees_c cluster[numb_cluster]; // stack problem related (limited)

	genrdm(cluster, numb_cluster);

	cout <<endl;
	for(unsigned i=0; i<numb_cluster; i++){
		cout << "Cluster nbr : " << i << " -> x = " << cluster[i].x << " y = " << cluster[i].y << endl;
	}
	cout <<endl;

	const char *scriptDirectoryName = "/home/rishie/Documents/Project/src/python/";

 	Py_Initialize();

    PyObject *sysPath = PySys_GetObject((char *)"path");
    PyObject *path = PyString_FromString(scriptDirectoryName);
    int result = PyList_Insert(sysPath, 0, path);
    pModule = PyImport_ImportModule("plot_xy");

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv[1]);
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(argc - 2);
            for (i = 0; i < argc - 2; ++i) {
                pValue = PyInt_FromLong(atoi(argv[i + 2]));
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                /* pValue reference stolen here: */
                PyTuple_SetItem(pArgs, i, pValue);
            }
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pValue != NULL) {
                printf("Result of call: %ld\n", PyInt_AsLong(pValue));
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", argv[1]);
        return 1;
    }
    Py_Finalize();



	cout << "The end!" << endl;
	return 0;
}