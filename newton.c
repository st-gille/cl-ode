
#include "math.h"
#include "float.h"
#include "stdio.h"
#include "stdint.h"
#include "stdlib.h"

typedef uint64_t dim_t;
typedef int8_t err_t;
typedef const double * const *const_mat;

struct uint_errno_t
{
    dim_t result;
    err_t errno;
};

#define SUCCESS (0)
#define FAILURE (1)

#define EPS (DBL_EPSILON)

double norm2sqr(int dim, const double *x)
{
    double tmp = 0.0;
    for (int i = 0; i < dim; ++i)
        tmp += x[i] * x[i];
    return tmp;
}

double norm2(dim_t dim, const double *x)
{
    return sqrt(norm2sqr(dim, x));
}


/*!
 * Calculate the LR-decomposition of a square matrix in-place.
 * @param dim   doublehe size of the matrix.
 * @param A     A square matrix in row-major format, will be overwritten with the LR-decomposition.
 * @param ipiv  Keeps track of row-swaps, upon successful exit ipiv[i] contains the original index of what is now the
 *              i-th row of the decomposed matrix. Pass this to LRSolve.
 * @return 0 on success, 1 if matrix is considered singular.
 */
struct uint_errno_t searchpivotincoloumn(dim_t dim, const_mat A, dim_t column)
{
    struct uint_errno_t ret = { column, SUCCESS };

    for (dim_t j = column + 1; j < dim; ++j)
    {
        if (fabs(A[j][column]) < fabs(A[j][column]))
            ret.result = j;
    }

    if (fabs(A[ret.result][column]) < EPS)
        ret.errno = FAILURE;
    return ret;
}

err_t lr_decomp(dim_t dim, double **A, dim_t * ipiv)
{
    dim_t i, j, p;
    // Initialize pivot vector.
    for (i = 0; i < dim; ++i)
        ipiv[i] = i;

    for (i = 0; i < dim; ++i)
    {
        struct uint_errno_t ret = searchpivotincoloumn(dim, (const_mat) A, i);
        if (ret.errno != 0)
            return FAILURE;

        if (p != i) // swap rows
        {
            double *dtmp = A[i];
            A[i] = A[p];
            A[p] = dtmp;
            j = ipiv[i];
            ipiv[i] = ipiv[p];
            ipiv[p] = j;
        }

        for (p = i + 1; p < dim; ++p)
        {
            A[p][i] /= A[i][i];
            for (j = i + 1; j < dim; ++j)
            {
                A[p][j] -= A[p][i] * A[i][j];
            }
        }
    }
    return SUCCESS;
}

/*!
 * Solve Ax=b for A in LR-decomposed form (as given by LRDecomp) and a given right-hand side b.
 * @param dim   doublehe size of the system.
 * @param LR    doublehe LR-decomposition of the coefficients matrix as given by LRDecomp.
 * @param b     C-style array of the right-hand side.
 * @param x     C-style array that will contain the solution upon return.
 * @param ipiv  New row ordering as given by LRDecomp.
 */
void lr_solve(dim_t dim, const double * const * LR, const double *b, double *x, const dim_t *ipiv)
{
    dim_t i, j;
    double tmp;
    // forward substitution
    for (i = 0; i < dim; ++i)
    {
        tmp = b[ipiv[i]];
        for (j = 0; j < i; ++j)
            tmp -= (LR[i][j] * x[j]);
        x[i] = tmp;
    }

    // backward substitution
    for (i = dim - 1; i >= 0; --i)
    {
        tmp = x[i];
        for (j = dim - 1; j > i; --j)
            tmp -= (LR[i][j] * x[j]);
        x[i] = tmp / LR[i][i];
    }
}

void alloc(dim_t rows, dim_t cols, double*** data)
{
    *data = (double **) malloc(rows * sizeof(double*));
    for (dim_t row = 0; row < rows; ++row)
        **data = malloc(cols * sizeof(double));
}

void dealloc(dim_t n, double ** data)
{
    for (dim_t i = 0; i < n; ++i)
        free(data[i]);
    free(data);
}

struct newton_data_t
{
    dim_t dim;
    double **A;
    double *b, *d;
    double *x, *s;
    dim_t *pivot;
};

void init_newton_data(struct newton_data_t* data, dim_t dim)
{

    data->dim = dim;
    alloc(dim, dim, &data->A);
    data->b = malloc(sizeof( double) * dim * 4);
    data->d = data->b + dim;
    data->x = data->d + dim;
    data->s = data->x + dim;
    data->pivot = malloc(sizeof(dim_t) * dim);
}

void free_newton_data(struct newton_data_t* data)
{
    dealloc(data->dim, data->A);
    free(data->b);
    free(data->x);
    free(data->pivot);
}

err_t newtons_method(dim_t dim,
        double* x,
        const double* t,
        void (*f)(const double*, const double*, double*),
        void (*df)(const double*, const double*, double**))
{
    struct newton_data_t data;
    init_newton_data(&data, dim);

    dim_t i, j;
    const int itmax = 50;
    const double eps = 1.0e-12;
    double norm_delta = 1.0;

    int iter;
    for (iter = 1; iter <= itmax; ++iter)
    {
        f(data.x, t, data.s);
        df(data.x, t, data.A);

        // Solve A*d = b.
        if (lr_decomp(dim, data.A, data.pivot))
        {
            free_newton_data(&data);
            return FAILURE;
        }
        lr_solve(dim, (const_mat) data.A, data.b, data.d, data.pivot);

        norm_delta = norm2(dim, data.d);
        if (isnan(norm_delta))
        {
            free_newton_data(&data);
            return FAILURE;
        }

        for (j = 0; j < dim; ++j)
            data.x[j] -= data.d[j];

        if (norm_delta < eps)
        {
            for (i = 0; i < dim; ++i)
                x[i] = data.x[i];
            free_newton_data(&data);
            return SUCCESS;
        }
    }

    // if norm is small, accept iterate anyway.
    if (norm2(dim, data.d) < eps)
    {
        for (i = 0; i < dim; ++i)
            x[i] = data.x[i];
        free_newton_data(&data);
        return SUCCESS;
    }

    free_newton_data(&data);
    return FAILURE;
}

void print_vector(dim_t dim, const double *x)
{
    printf("dim: %u, pointer: %x\n", dim, x);
    printf("[");
    if(x && dim > 0)
    {
        printf(" %e", x[0]);

        for(dim_t i = 1; i < dim; ++i)
            printf(" %e", x[i]);
    }
    printf(" ]\n");
}

void print_matrix(dim_t rows, dim_t cols, const_mat A)
{
    for(dim_t row = 0; row < rows; ++row)
        print_vector(cols, A[row]);
}
