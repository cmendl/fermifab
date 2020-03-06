#include "tensor_op.h"
#include "util.h"
#include "binio.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Uniform distance (infinity norm) between 'x' and 'y'
///
static inline double UniformDistance(const int n, const double *x, const double *y)
{
	double d = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		d = fmax(d, fabs(x[i] - y[i]));
	}

	return d;
}


//________________________________________________________________________________________________________________________
//


int TestTensorOp(void)
{
	printf("Testing tensor_op module...\n");

	int status;
	double err = 0;

	const int orbs = 6;

	double *A = (double *)malloc(orbs*orbs * sizeof(double));
	status = ReadData("../test/test_tensor_op_A.dat", A, sizeof(double), orbs*orbs);
	if (status < 0) { return status; }

	// determinant
	{
		// copy 'A' since determinant function overwrites it
		double *Acopy = (double *)malloc(orbs*orbs * sizeof(double));
		memcpy(Acopy, A, orbs*orbs * sizeof(double));

		double d = Det(orbs, Acopy);

		free(Acopy);

		// load reference determinant from disk
		double d_ref = 0;
		status = ReadData("../test/test_tensor_op_detA.dat", &d_ref, sizeof(double), 1);
		if (status < 0) { return status; }

		err += fabs(d - d_ref) / fabs(d_ref);
	}

	// tensor product
	{
		const int N = 4;

		sparse_array_t AN = { 0 };
		status = TensorOp(orbs, N, A, &AN);
		if (status < 0) { return status; }

		const int nelem = IntegerProduct(AN.dims, AN.rank);

		double *ANd = (double *)malloc(nelem * sizeof(double));
		if (ANd == NULL) { return -1; }
		SparseToDense(&AN, ANd);

		// load reference matrix from disk
		double *AN_ref = (double *)malloc(nelem * sizeof(double));
		if (AN_ref == NULL) { return -1; }
		status = ReadData("../test/test_tensor_op_AN.dat", AN_ref, sizeof(double), nelem);
		if (status < 0) { return status; }

		err += UniformDistance(nelem, ANd, AN_ref);

		free(AN_ref);
		free(ANd);
		DeleteSparseArray(&AN);
	}

	// clean up
	free(A);

	return (err < 1e-12 ? 0 : 1);
}
