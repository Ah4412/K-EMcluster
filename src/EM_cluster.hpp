#ifndef EM_cluster_hpp
#define EM_cluster_hpp

void EM_cluster_p(int nbr_cluster, int size)
{	
	PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue, *pXVec;

    const char *scriptDirectoryName = "/home/alex/Documents/Project/src/python/";

 	Py_Initialize();

    PyObject *sysPath = PySys_GetObject((char *)"path");
    PyObject *path = PyString_FromString(scriptDirectoryName);
    int result = PyList_Insert(sysPath, 0, path);
    pModule = PyImport_ImportModule("rdm_pnt");

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, "rdm_pnt");
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {

            pArgs = PyTuple_New(1);

            pXVec = PyTuple_New(2);

                pValue = PyInt_FromLong((long)nbr_cluster);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument x\n");
                    return;
                }
                PyTuple_SetItem(pXVec, 0, pValue);
                //set Y's value on other array.
                pValue = PyInt_FromLong((long)size);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument y\n");
                    return;
                }
                PyTuple_SetItem(pXVec, 1, pValue);  

            PyTuple_SetItem(pArgs, 0, pXVec);
			
            pValue = PyObject_CallObject(pFunc, pArgs);
            //Py_INCREF(pArgs);
            Py_DECREF(pArgs);

            if (PyInt_AsLong(pValue) == 1) {
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
    //Py_Finalize();
}

#endif