#ifndef plot_dyn_hpp
#define plot_dyn_hpp

void plot_dyn_p(const vector<coordonnees_c> &c, const coordonnees_p *p, const int numb_cluster)
{	
	PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue, *pXVec, *pYVec;
    PyObject *pX, *pY;
    int i;



    const char *scriptDirectoryName = "/home/alex/Documents/Project/src/python/";

 	Py_Initialize();

    PyObject *sysPath = PySys_GetObject((char *)"path");
    PyObject *path = PyString_FromString(scriptDirectoryName);
    int result = PyList_Insert(sysPath, 0, path);
    pModule = PyImport_ImportModule("plot");

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, "plot");
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {

            pArgs = PyTuple_New(4);

            pXVec = PyTuple_New(numb_cluster); 
            pYVec = PyTuple_New(numb_cluster);
            pX = PyTuple_New(table_size);
            pY = PyTuple_New(table_size);

            for (i = 0; i < numb_cluster; i++){
                pValue = PyFloat_FromDouble(c[i].x);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument x\n");
                    return;
                }
                PyTuple_SetItem(pXVec, i, pValue);
                //set Y's value on other array.
                pValue = PyFloat_FromDouble(c[i].y);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument y\n");
                    return;
                }
                PyTuple_SetItem(pYVec, i, pValue);   
            }

            for (i = 0; i < table_size; i++) {
                pValue = PyFloat_FromDouble(p[i].x);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument x\n");
                    return;
                }
                PyTuple_SetItem(pX, i, pValue);
                //set Y's value on other array.
                pValue = PyFloat_FromDouble(p[i].y);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument y\n");
                    return;
                }
                PyTuple_SetItem(pY, i, pValue);   
            }

            PyTuple_SetItem(pArgs, 0, pXVec);
			PyTuple_SetItem(pArgs, 1, pYVec);
            PyTuple_SetItem(pArgs, 2, pX);
            PyTuple_SetItem(pArgs, 3, pY);
			
            pValue = PyObject_CallObject(pFunc, pArgs);
            //Py_INCREF(pArgs);
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
                return;
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
        return;
    }
    Py_Finalize();
}

#endif