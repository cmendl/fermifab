#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "fermi_map.h"
#include "generate_rdm.h"
#include "tensor_op.h"
#include <stdbool.h>
#include <inttypes.h>


//________________________________________________________________________________________________________________________
//


static PyObject *fermi2coords(PyObject *self, PyObject *args)
{
	// suppress "unused parameter" warning
	(void)self;

	PyObject *obj_orbs;     // list of number of orbitals
	PyObject *obj_N;        // particle numbers

	if (!PyArg_ParseTuple(args, "OO", &obj_orbs, &obj_N)) {
		PyErr_SetString(PyExc_SyntaxError, "error parsing input; syntax: fermi2coords(orbs, N)");
		return NULL;
	}

	int _orbs[64] = { 0 };
	int _N[64]    = { 0 };
	fermi_config_t config = {
		.orbs = _orbs,
		.N = _N,
		.nc = 0
	};

	// 'orbs'
	{
		PyObject *iterator = PyObject_GetIter(obj_orbs);
		if (iterator == NULL)
		{
			PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot iterate over 'orbs' argument; syntax: fermi2coords(orbs, N)");
			return NULL;
		}

		PyObject *item;
		while ((item = PyIter_Next(iterator)))
		{
			long n = PyLong_AsLong(item);
			if (PyErr_Occurred()) {
				PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot interpret 'orbs' item as integer; syntax: fermi2coords(orbs, N)");
				Py_DECREF(item);
				Py_DECREF(iterator);
				return NULL;
			}
			Py_DECREF(item);

			if (n <= 0) {
				PyErr_SetString(PyExc_ValueError, "entries in 'orbs' argument must be positive; syntax: fermi2coords(orbs, N)");
				Py_DECREF(iterator);
				return NULL;
			}
			else if (config.nc == 64) {
				PyErr_SetString(PyExc_ValueError, "too many entries in 'orbs' argument; syntax: fermi2coords(orbs, N)");
				Py_DECREF(iterator);
				return NULL;
			}

			config.orbs[config.nc] = (int)n;
			config.nc++;
		}

		Py_DECREF(iterator);
	}

	// 'N'
	{
		int count = 0;

		PyObject *iterator = PyObject_GetIter(obj_N);
		if (iterator == NULL)
		{
			PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot iterate over 'N' argument; syntax: fermi2coords(orbs, N)");
			return NULL;
		}

		PyObject *item;
		while ((item = PyIter_Next(iterator)))
		{
			long n = PyLong_AsLong(item);
			if (PyErr_Occurred()) {
				PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot interpret 'N' item as integer; syntax: fermi2coords(orbs, N)");
				Py_DECREF(item);
				Py_DECREF(iterator);
				return NULL;
			}
			Py_DECREF(item);

			if (n < 0) {
				PyErr_SetString(PyExc_ValueError, "entries in 'N' argument must be non-negative; syntax: fermi2coords(orbs, N)");
				Py_DECREF(iterator);
				return NULL;
			}
			else if (count == 64) {
				PyErr_SetString(PyExc_ValueError, "too many entries in 'N' argument; syntax: fermi2coords(orbs, N)");
				Py_DECREF(iterator);
				return NULL;
			}

			config.N[count] = (int)n;
			count++;
		}

		Py_DECREF(iterator);

		if (count != config.nc) {
			PyErr_SetString(PyExc_SyntaxError, "number of items in 'orbs' and 'N' must be the same; syntax: fermi2coords(orbs, N)");
			return NULL;
		}
	}
	int Ntot = IntegerSum(config.N, config.nc);

	int i;
	for (i = 0; i < config.nc; i++)
	{
		if (config.N[i] > config.orbs[i]) {
			PyErr_SetString(PyExc_ValueError, "particle number cannot be larger than the corresponding number of orbitals; syntax: fermi2coords(orbs, N)");
			return NULL;
		}
	}

	fermi_map_t baseMap;
	int status = FermiMap(&config, &baseMap);
	if (status < 0) {
		PyErr_SetString(PyExc_RuntimeError, "internal error occurred, probably out of memory");
		return NULL;
	}
	fermi_coords_t *x = (fermi_coords_t *)malloc(Ntot*sizeof(fermi_coords_t));

	npy_intp dims[2] = { baseMap.num, Ntot };
	PyArrayObject *coords_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_INT);
	if (coords_arr == NULL) {
		PyErr_SetString(PyExc_RuntimeError, "error creating to-be-returned matrix of coordinates");
		return NULL;
	}

	int *coords = (int *)PyArray_DATA(coords_arr);
	for (i = 0; i < baseMap.num; i++)
	{
		// get coordinates
		FermiDecode(baseMap.map[i], x, Ntot);

		int j;
		for (j = 0; j < Ntot; j++)
		{
			coords[i*Ntot + j] = x[j];
		}
	}

	// clean up
	free(x);
	free(baseMap.map);

	return (PyObject *)coords_arr;
}


//________________________________________________________________________________________________________________________
//


static PyObject *gen_rdm(PyObject *self, PyObject *args)
{
	// suppress "unused parameter" warning
	(void)self;

	PyObject *obj_orbs;     // list of number of orbitals
	PyObject *obj_p1;       // 'p1' particle numbers
	PyObject *obj_N1;       // 'N1' particle numbers
	PyObject *obj_N2;       // 'N1' particle numbers

	if (!PyArg_ParseTuple(args, "OOOO", &obj_orbs, &obj_p1, &obj_N1, &obj_N2)) {
		PyErr_SetString(PyExc_SyntaxError, "error parsing input; syntax: gen_rdm(orbs, p1, N1, N2)");
		return NULL;
	}

	int orbs[64] = { 0 };
	int p1[64]   = { 0 };
	int N1[64]   = { 0 };
	int N2[64]   = { 0 };
	int nc = 0;

	// 'orbs'
	{
		PyObject *iterator = PyObject_GetIter(obj_orbs);
		if (iterator == NULL)
		{
			PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot iterate over 'orbs' argument; syntax: gen_rdm(orbs, p1, N1, N2)");
			return NULL;
		}

		PyObject *item;
		while ((item = PyIter_Next(iterator)))
		{
			long n = PyLong_AsLong(item);
			if (PyErr_Occurred()) {
				PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot interpret 'orbs' item as integer; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(item);
				Py_DECREF(iterator);
				return NULL;
			}
			Py_DECREF(item);

			if (n < 0) {
				PyErr_SetString(PyExc_ValueError, "entries in 'orbs' argument must be non-negative; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(iterator);
				return NULL;
			}
			else if (nc == 64) {
				PyErr_SetString(PyExc_ValueError, "too many entries in 'orbs' argument; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(iterator);
				return NULL;
			}

			orbs[nc] = (int)n;
			nc++;
		}

		Py_DECREF(iterator);
	}

	// 'p1'
	{
		int count = 0;

		PyObject *iterator = PyObject_GetIter(obj_p1);
		if (iterator == NULL)
		{
			PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot iterate over 'p1' argument; syntax: gen_rdm(orbs, p1, N1, N2)");
			return NULL;
		}

		PyObject *item;
		while ((item = PyIter_Next(iterator)))
		{
			long n = PyLong_AsLong(item);
			if (PyErr_Occurred()) {
				PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot interpret 'p1' item as integer; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(item);
				Py_DECREF(iterator);
				return NULL;
			}
			Py_DECREF(item);

			if (n < 0) {
				PyErr_SetString(PyExc_ValueError, "entries in 'p1' argument must be non-negative; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(iterator);
				return NULL;
			}
			else if (count == 64) {
				PyErr_SetString(PyExc_ValueError, "too many entries in 'p1' argument; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(iterator);
				return NULL;
			}

			p1[count] = (int)n;
			count++;
		}

		Py_DECREF(iterator);

		if (count != nc) {
			PyErr_SetString(PyExc_SyntaxError, "number of items in 'orbs' and 'p1' must be the same; syntax: gen_rdm(orbs, p1, N1, N2)");
			return NULL;
		}
	}

	// 'N1'
	{
		int count = 0;

		PyObject *iterator = PyObject_GetIter(obj_N1);
		if (iterator == NULL)
		{
			PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot iterate over 'N1' argument; syntax: gen_rdm(orbs, p1, N1, N2)");
			return NULL;
		}

		PyObject *item;
		while ((item = PyIter_Next(iterator)))
		{
			long n = PyLong_AsLong(item);
			if (PyErr_Occurred()) {
				PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot interpret 'N1' item as integer; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(item);
				Py_DECREF(iterator);
				return NULL;
			}
			Py_DECREF(item);

			if (n < 0) {
				PyErr_SetString(PyExc_ValueError, "entries in 'N1' argument must be non-negative; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(iterator);
				return NULL;
			}
			else if (count == 64) {
				PyErr_SetString(PyExc_ValueError, "too many entries in 'N1' argument; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(iterator);
				return NULL;
			}

			N1[count] = (int)n;
			count++;
		}

		Py_DECREF(iterator);

		if (count != nc) {
			PyErr_SetString(PyExc_SyntaxError, "number of items in 'orbs' and 'N1' must be the same; syntax: gen_rdm(orbs, p1, N1, N2)");
			return NULL;
		}
	}

	// 'N2'
	{
		int count = 0;

		PyObject *iterator = PyObject_GetIter(obj_N2);
		if (iterator == NULL)
		{
			PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot iterate over 'N2' argument; syntax: gen_rdm(orbs, p1, N1, N2)");
			return NULL;
		}

		PyObject *item;
		while ((item = PyIter_Next(iterator)))
		{
			long n = PyLong_AsLong(item);
			if (PyErr_Occurred()) {
				PyErr_SetString(PyExc_SyntaxError, "error parsing input: cannot interpret 'N2' item as integer; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(item);
				Py_DECREF(iterator);
				return NULL;
			}
			Py_DECREF(item);

			if (n < 0) {
				PyErr_SetString(PyExc_ValueError, "entries in 'N2' argument must be non-negative; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(iterator);
				return NULL;
			}
			else if (count == 64) {
				PyErr_SetString(PyExc_ValueError, "too many entries in 'N2' argument; syntax: gen_rdm(orbs, p1, N1, N2)");
				Py_DECREF(iterator);
				return NULL;
			}

			N2[count] = (int)n;
			count++;
		}

		Py_DECREF(iterator);

		if (count != nc) {
			PyErr_SetString(PyExc_SyntaxError, "number of items in 'orbs' and 'N2' must be the same; syntax: gen_rdm(orbs, p1, N1, N2)");
			return NULL;
		}
	}

	int i;
	for (i = 0; i < nc; i++)
	{
		if (N1[i] - N2[i] + p1[i] < 0) {
			PyErr_SetString(PyExc_ValueError, "all entries in 'N1 - N2 + p1' must be non-negative; syntax: gen_rdm(orbs, p1, N1, N2)");
			return NULL;
		}
	}

	// actually compute kernel tensor
	sparse_array_t K = { 0 };
	int status = GenerateRDM(orbs, p1, N1, N2, nc, &K);
	if (status < 0) {
		PyErr_SetString(PyExc_RuntimeError, "internal error occurred, probably out of memory");
		DeleteSparseArray(&K);
		return NULL;
	}

	// kernel tensor dimensions
	assert(K.rank == 4);
	PyObject *dims_obj = Py_BuildValue("(iiii)", K.dims[0], K.dims[1], K.dims[2], K.dims[3]);

	// construct array of values
	npy_intp dims_val[1] = { K.nnz };
	PyArrayObject *val_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims_val, NPY_DOUBLE);
	if (val_arr == NULL) {
		PyErr_SetString(PyExc_RuntimeError, "error creating to-be-returned value vector");
		Py_DECREF(dims_obj);
		DeleteSparseArray(&K);
		return NULL;
	}
	memcpy(PyArray_DATA(val_arr), K.val, K.nnz * sizeof(double));

	// construct array of indices
	npy_intp dims_ind[2] = { K.nnz, K.rank };
	PyArrayObject *ind_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims_ind, sizeof(K.ind[0]) == 4 ? NPY_INT32 : NPY_INT64);
	if (ind_arr == NULL) {
		PyErr_SetString(PyExc_RuntimeError, "error creating to-be-returned array of indices");
		Py_DECREF(dims_obj);
		Py_DECREF(val_arr);
		DeleteSparseArray(&K);
		return NULL;
	}
	memcpy(PyArray_DATA(ind_arr), K.ind, K.nnz*K.rank * sizeof(K.ind[0]));

	// clean up
	DeleteSparseArray(&K);

	return Py_BuildValue("(OOO)", dims_obj, val_arr, ind_arr);
}


//________________________________________________________________________________________________________________________
//


static PyObject *tensor_op(PyObject *self, PyObject *args)
{
	// suppress "unused parameter" warning
	(void)self;

	PyObject *Ain;
	int N;
	if (!PyArg_ParseTuple(args, "Oi", &Ain, &N)) {
		PyErr_SetString(PyExc_SyntaxError, "error parsing input; syntax: tensor_op(A, N)");
		return NULL;
	}

	if (N <= 0) {
		PyErr_SetString(PyExc_ValueError, "'N' must be positive; syntax: tensor_op(A, N)");
		return NULL;
	}

	// find out if we should aim for a real or complex matrix
	bool use_complex;
	{
		PyArrayObject *arr = NULL;
		PyArray_Descr *dtype = NULL;
		int ndim = 0;
		npy_intp dims[NPY_MAXDIMS];
		if (PyArray_GetArrayParamsFromObject(Ain, /*requested_dtype*/NULL, /*writeable*/false, &dtype, &ndim, dims, &arr, NULL) < 0)
		{
			PyErr_SetString(PyExc_SyntaxError, "cannot get array parameters of 'A'; syntax: tensor_op(A, N)");
			return NULL;
		}
		if (arr == NULL)
		{
			PyErr_SetString(PyExc_SyntaxError, "cannot interpret 'A' as array; syntax: tensor_op(A, N)");
			return NULL;
		}

		PyArray_Descr *desc = PyArray_DESCR(arr);
		use_complex = (desc->kind == 'c');

		Py_DECREF(arr);
	}

	if (!use_complex)
	{
		PyArrayObject *A = (PyArrayObject *)PyArray_ContiguousFromObject(Ain, NPY_DOUBLE, 2, 2);
		if (A == NULL)
		{
			PyErr_SetString(PyExc_ValueError, "cannot interpret 'A' as real matrix");
			return NULL;
		}

		if (PyArray_DIM(A, 0) != PyArray_DIM(A, 1))
		{
			PyErr_SetString(PyExc_ValueError, "'A' must be a square matrix");
			Py_DECREF(A);
			return NULL;
		}

		const int orbs = PyArray_DIM(A, 0);

		if (N > orbs) {
			PyErr_SetString(PyExc_ValueError, "'N' cannot be larger than number of orbitals; syntax: tensor_op(A, N)");
			Py_DECREF(A);
			return NULL;
		}

		sparse_array_t AN = { 0 };
		int status = TensorOp(orbs, N, PyArray_DATA(A), &AN);
		if (status < 0) {
			PyErr_SetString(PyExc_RuntimeError, "internal error occurred, probably out of memory");
			DeleteSparseArray(&AN);
			Py_DECREF(A);
			return NULL;
		}

		Py_DECREF(A);

		// dimensions
		assert(AN.rank == 2);
		PyObject *dims_obj = Py_BuildValue("(ii)", AN.dims[0], AN.dims[1]);

		// construct array of values
		npy_intp dims_val[1] = { AN.nnz };
		PyArrayObject *val_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims_val, NPY_DOUBLE);
		if (val_arr == NULL) {
			PyErr_SetString(PyExc_RuntimeError, "error creating to-be-returned value vector");
			Py_DECREF(dims_obj);
			DeleteSparseArray(&AN);
			return NULL;
		}
		memcpy(PyArray_DATA(val_arr), AN.val, AN.nnz * sizeof(double));

		// construct array of indices
		npy_intp dims_ind[2] = { AN.nnz, AN.rank };
		PyArrayObject *ind_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims_ind, sizeof(AN.ind[0]) == 4 ? NPY_INT32 : NPY_INT64);
		if (ind_arr == NULL) {
			PyErr_SetString(PyExc_RuntimeError, "error creating to-be-returned array of indices");
			Py_DECREF(dims_obj);
			Py_DECREF(val_arr);
			DeleteSparseArray(&AN);
			return NULL;
		}
		memcpy(PyArray_DATA(ind_arr), AN.ind, AN.nnz*AN.rank * sizeof(AN.ind[0]));

		// clean up
		DeleteSparseArray(&AN);

		return Py_BuildValue("(OOO)", dims_obj, val_arr, ind_arr);
	}
	else // complex-valued input matrix
	{
		PyArrayObject *A = (PyArrayObject *)PyArray_ContiguousFromObject(Ain, NPY_CDOUBLE, 2, 2);
		if (A == NULL)
		{
			PyErr_SetString(PyExc_ValueError, "cannot interpret 'A' as complex matrix");
			return NULL;
		}

		if (PyArray_DIM(A, 0) != PyArray_DIM(A, 1))
		{
			PyErr_SetString(PyExc_ValueError, "'A' must be a square matrix");
			Py_DECREF(A);
			return NULL;
		}

		const int orbs = PyArray_DIM(A, 0);

		if (N > orbs) {
			PyErr_SetString(PyExc_ValueError, "'N' cannot be larger than number of orbitals; syntax: tensor_op(A, N)");
			Py_DECREF(A);
			return NULL;
		}

		sparse_complex_array_t AN = { 0 };
		int status = TensorOpComplex(orbs, N, PyArray_DATA(A), &AN);
		if (status < 0) {
			PyErr_SetString(PyExc_RuntimeError, "internal error occurred, probably out of memory");
			DeleteSparseComplexArray(&AN);
			Py_DECREF(A);
			return NULL;
		}

		Py_DECREF(A);

		// dimensions
		assert(AN.rank == 2);
		PyObject *dims_obj = Py_BuildValue("(ii)", AN.dims[0], AN.dims[1]);

		// construct array of values
		npy_intp dims_val[1] = { AN.nnz };
		PyArrayObject *val_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims_val, NPY_CDOUBLE);
		if (val_arr == NULL) {
			PyErr_SetString(PyExc_RuntimeError, "error creating to-be-returned value vector");
			Py_DECREF(dims_obj);
			DeleteSparseComplexArray(&AN);
			return NULL;
		}
		memcpy(PyArray_DATA(val_arr), AN.val, AN.nnz * sizeof(double complex));

		// construct array of indices
		npy_intp dims_ind[2] = { AN.nnz, AN.rank };
		PyArrayObject *ind_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims_ind, sizeof(AN.ind[0]) == 4 ? NPY_INT32 : NPY_INT64);
		if (ind_arr == NULL) {
			PyErr_SetString(PyExc_RuntimeError, "error creating to-be-returned array of indices");
			Py_DECREF(dims_obj);
			Py_DECREF(val_arr);
			DeleteSparseComplexArray(&AN);
			return NULL;
		}
		memcpy(PyArray_DATA(ind_arr), AN.ind, AN.nnz*AN.rank * sizeof(AN.ind[0]));

		// clean up
		DeleteSparseComplexArray(&AN);

		return Py_BuildValue("(OOO)", dims_obj, val_arr, ind_arr);
	}
}


//________________________________________________________________________________________________________________________
//


static PyMethodDef methods[] = {
	{ "fermi2coords", fermi2coords, METH_VARARGS, "Enumerate all N-particle Slater basis states for 'orbs' available orbitals." },
	{ "gen_rdm",      gen_rdm,      METH_VARARGS, "Generate sparse kernel tensor for computing reduced density matrices." },
	{ "tensor_op",    tensor_op,    METH_VARARGS, "Matrix representation of the N-fold tensor product of an operator." },
	{ NULL, NULL, 0, NULL }     // sentinel
};


static struct PyModuleDef module = {
	PyModuleDef_HEAD_INIT,
	"fermifab.kernel",      // name of module
	NULL,                   // module documentation, may be NULL
	-1,                     // size of per-interpreter state of the module, or -1 if the module keeps state in global variables
	methods                 // module methods
};


PyMODINIT_FUNC PyInit_kernel(void)
{
	// import NumPy array module (required)
	import_array();

	PyObject *m = PyModule_Create(&module);
	if (m == NULL) {
		return NULL;
	}

	return m;
}
