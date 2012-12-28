#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "newton.h"

#ifdef NDEBUG
#   define DEBUG_PRINT(...)
#   define DEBUG(X)
#else
#   define DEBUG_PRINT printf
#   define DEBUG(X) X
#endif

#define DEFINE_PRINT_VECTOR_T(type, suffix, format)             \
    void print_vector##suffix(dim_t dim, const type *x)         \
{                                                           \
    printf("[");                                            \
    if(x && dim > 0)                                        \
    {                                                       \
        for(dim_t i = 0; i < dim; ++i)                      \
        printf(format, x[i]);                               \
    }                                                       \
    printf(" ]\n");                                         \
}

DEFINE_PRINT_VECTOR_T(int, _i, " %d")
DEFINE_PRINT_VECTOR_T(dim_t, _u, " %u")
DEFINE_PRINT_VECTOR_T(double, , " %e")
DEFINE_PRINT_VECTOR_T(double *, _x, " %x")

#undef DEFINE_PRINT_VECTOR_T

void print_matrix(dim_t rows, dim_t cols, const_mat A)
{
    for(dim_t row = 0; row < rows; ++row)
        print_vector(cols, A[row]);
}

double norm2sqr(dim_t dim, const double *x)
{
    double tmp = 0.0;
    for (dim_t i = 0; i < dim; ++i)
        tmp += x[i] * x[i];
    DEBUG_PRINT("norm2sqr: dim = %u, res = %f, x =", dim, tmp);

    return tmp;
}

double norm2(dim_t dim, const double *x)
{
    DEBUG_PRINT("in norm2: dim = %u, x = ", dim);
    DEBUG(print_vector(dim, x));
    return sqrt(norm2sqr(dim, x));
}

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
    dim_t i, j;
    double dtmp;
    struct uint_errno_t ret;

    for (i = 0; i < dim; ++i)
        ipiv[i] = i;

    for (i = 0; i < dim; ++i)
    {
        ret = searchpivotincoloumn(dim, (const_mat) A, i);
        if (ret.errno == FAILURE)
            return FAILURE;

        if (ret.result != i)
        {
            for(j = 0; j < dim; ++j)
            {
                dtmp = A[i][j];
                A[i][j] = A[ret.result][j];
                A[ret.result][j] = dtmp;
            }
            j = ipiv[i];
            ipiv[i] = ipiv[ret.result];
            ipiv[ret.result] = j;
        }

        for (ret.result = i + 1; ret.result < dim; ++ret.result)
        {
            A[ret.result][i] /= A[i][i];
            for (j = i + 1; j < dim; ++j)
                A[ret.result][j] -= A[ret.result][i] * A[i][j];
        }
    }
    return SUCCESS;
}

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
    i = dim;
    do
    {
        --i;
        tmp = x[i];
        for (j = dim - 1; j > i; --j)
            tmp -= (LR[i][j] * x[j]);
        x[i] = tmp / LR[i][i];

    } while(i != 0);
}

void *matrix_alloc(dim_t rows, dim_t cols, dim_t size)
{
    DEBUG_PRINT("in alloc: %u %u %#x\n", rows, cols);
    void **data = malloc(rows * sizeof(void*));
    void *tmp = malloc(rows * cols * size);
    for (dim_t row = 0; row < rows; ++row)
    {
        data[row] = tmp;
        tmp += cols * size;
    }
    return data;
}

void matrix_dealloc(dim_t n, void **data)
{
    free(*data);
    free(data);
}

struct newton_data_t
{
    dim_t dim;
    double **A, *d, *s;
    dim_t *pivot;
};

void alloc_newton_data(struct newton_data_t* data, dim_t dim)
{
    data->dim = dim;
    DEBUG_PRINT("in alloc_newton_data\n");
    DEBUG_PRINT("data: %#x\n", data);
    data->A = matrix_alloc(dim, dim, sizeof(double));
    data->d = malloc(sizeof(double) * dim * 2);
    data->s = data->d + dim;
    data->pivot = malloc(sizeof(dim_t) * dim);
}

void free_newton_data(struct newton_data_t* data)
{
    matrix_dealloc(data->dim, (void **) data->A);
    free(data->d);
    free(data->pivot);
}

err_t newtons_method(dim_t dim,
        double* x,
        void (*f)(const double*, double*),
        void (*df)(const double*, double**))
{
    DEBUG_PRINT("in newtons_method: (dim x f df) = (%d %#x %#X %#X)\n", dim, x, f, df);
    struct newton_data_t data;
    alloc_newton_data(&data, dim);
    DEBUG_PRINT("allocated newton data: %#x\n", &data);

    dim_t j;
    const int itmax = 50;
    const double eps = 1.0e-12;
    double norm_delta = 1.0;

    DEBUG_PRINT("x = ");
    DEBUG(print_vector(dim, x));
    int iter;
    for (iter = 1; iter <= itmax; ++iter)
    {
        DEBUG_PRINT("----------- step %d----------------\nf(x) = ", iter);

        f(x, data.s);
        DEBUG(print_vector(dim, data.s));
        DEBUG_PRINT("|f(x)| = %e, ", norm2(dim, data.s));

        if(norm2(dim, data.s) < eps)
        {
            free_newton_data(&data);
            return SUCCESS;
        }

        df(x, data.A);
        DEBUG_PRINT("A = ");
        DEBUG(print_matrix(dim, dim, (const_mat) data.A));
        if (lr_decomp(dim, data.A, data.pivot))
        {
            free_newton_data(&data);
            return FAILURE;
        }
        lr_solve(dim, (const_mat) data.A, data.s, data.d, data.pivot);

        norm_delta = norm2(dim, data.d);
        DEBUG_PRINT("|delta| = %e, delta = ", norm_delta);
        DEBUG(print_vector(dim, data.d));

        if (isnan(norm_delta))
        {
            free_newton_data(&data);
            return FAILURE;
        }

        for (j = 0; j < dim; ++j)
            x[j] -= data.d[j];
        DEBUG_PRINT("x = ");
        DEBUG(print_vector(dim, x));

        if (norm_delta < eps)
        {
            free_newton_data(&data);
            return SUCCESS;
        }
    }

    free_newton_data(&data);
    return FAILURE;
}

err_t newtons_method_implicit(dim_t dim,
        double* x,
        const double* t,
        void (*f)(const double*, const double*, double*),
        void (*df)(const double*, const double*, double**))
{
    DEBUG_PRINT("in newtons_method: (dim x t f df) = (%d %#x %#x %#X %#X)\n", dim, x, t, f, df);
    struct newton_data_t data;
    alloc_newton_data(&data, dim);
    DEBUG_PRINT("allocated newton data: %#x\n", &data);

    dim_t j;
    const int itmax = 50;
    const double eps = 1.0e-12;
    double norm_delta = 1.0;

    DEBUG_PRINT("x = ");
    DEBUG(print_vector(dim, x));
    DEBUG_PRINT("t = ");
    DEBUG(print_vector(dim, t));
    int iter;
    for (iter = 1; iter <= itmax; ++iter)
    {
        DEBUG_PRINT("----------- step %d----------------\n", iter);

        f(x, t, data.s);
        DEBUG_PRINT("|f(x,t)| = %e, ", norm2(dim, data.s));
        DEBUG(print_vector(dim, data.s));

        if(norm2(dim, data.s) < eps)
        {
            free_newton_data(&data);
            return SUCCESS;
        }

        df(x, t, data.A);
        DEBUG_PRINT("A = ");
        DEBUG(print_matrix(dim, dim, (const_mat) data.A));
        if (lr_decomp(dim, data.A, data.pivot))
        {
            free_newton_data(&data);
            return FAILURE;
        }
        lr_solve(dim, (const_mat) data.A, data.s, data.d, data.pivot);

        norm_delta = norm2(dim, data.d);
        DEBUG_PRINT("|delta| = %e, delta = ", norm_delta);
        DEBUG(print_vector(dim, data.d));

        if (isnan(norm_delta))
        {
            free_newton_data(&data);
            return FAILURE;
        }

        for (j = 0; j < dim; ++j)
            x[j] -= data.d[j];
        DEBUG_PRINT("x = ");
        DEBUG(print_vector(dim, x));

        if (norm_delta < eps)
        {
            free_newton_data(&data);
            return SUCCESS;
        }
    }
    free_newton_data(&data);
    return FAILURE;
}

#undef DEBUG_PRINT
