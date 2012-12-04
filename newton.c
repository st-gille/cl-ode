#include "newton.h"

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
    for (i = dim - 1; ; --i)
    {
        tmp = x[i];
        for (j = dim - 1; j > i; --j)
            tmp -= (LR[i][j] * x[j]);
        x[i] = tmp / LR[i][i];

        if (i == 0)
            break;
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
    double **A, *d, *s;
    dim_t *pivot;
};

void alloc_newton_data(struct newton_data_t* data, dim_t dim)
{
    data->dim = dim;
    alloc(dim, dim, &(data->A));
    data->d = malloc(sizeof(double) * dim * 2);
    data->s = data->d + dim;
    data->pivot = malloc(sizeof(dim_t) * dim);
}

void free_newton_data(struct newton_data_t* data)
{
    dealloc(data->dim, data->A);
    free(data->d);
    free(data->pivot);
}

err_t newtons_method(dim_t dim,
        double* x,
        const double* t,
        void (*f)(const double*, const double*, double*),
        void (*df)(const double*, const double*, double**))
{
    struct newton_data_t data;
    alloc_newton_data(&data, dim);

    dim_t j;
    const int itmax = 50;
    const double eps = 1.0e-12;
    double norm_delta = 1.0;

    int iter;
    for (iter = 1; iter <= itmax; ++iter)
    {

        f(x, t, data.s);

        if(norm2(dim, data.s) < eps)
        {
            free_newton_data(&data);
            return SUCCESS;
        }

        df(x, t, data.A);
        if (lr_decomp(dim, data.A, data.pivot))
        {
            free_newton_data(&data);
            return FAILURE;
        }
        lr_solve(dim, (const_mat) data.A, data.s, data.d, data.pivot);

        norm_delta = norm2(dim, data.d);
        if (isnan(norm_delta))
        {
            free_newton_data(&data);
            return FAILURE;
        }

        for (j = 0; j < dim; ++j)
            x[j] -= data.d[j];

        if (norm_delta < eps)
        {
            free_newton_data(&data);
            return SUCCESS;
        }
    }

    free_newton_data(&data);
    return FAILURE;
}

#define DEFINE_PRINT_VECTOR_T(type, suffix, format)             \
    void print_vector##suffix(dim_t dim, const type *x)         \
    {                                                           \
        printf("[");                                            \
        if(x && dim > 0)                                        \
        {                                                       \
            printf(format, x[0]);                               \
            for(dim_t i = 1; i < dim; ++i)                      \
            printf(format, x[i]);                               \
        }                                                       \
        printf(" ]\n");                                         \
    }

DEFINE_PRINT_VECTOR_T(int, _i, " %d")
DEFINE_PRINT_VECTOR_T(dim_t, _u, " %u")
DEFINE_PRINT_VECTOR_T(double, , " %e")

#undef DEFINE_PRINT_VECTOR_T

void print_matrix(dim_t rows, dim_t cols, const_mat A)
{
    for(dim_t row = 0; row < rows; ++row)
        print_vector(cols, A[row]);
}
