#include <stdint.h>

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

double norm2sqr(dim_t dim, const double *x);
double norm2(dim_t dim, const double *x);


/*!
 * Calculate the LR-decomposition of a square matrix in-place.
 * @param dim   doublehe size of the matrix.
 * @param A     A square matrix in row-major format, will be overwritten with the LR-decomposition.
 * @param ipiv  Keeps track of row-swaps, upon successful exit ipiv[i] contains the original index of what is now the
 *              i-th row of the decomposed matrix. Pass this to LRSolve.
 * @return 0 on success, 1 if matrix is considered singular.
 */
err_t lr_decomp(dim_t dim, double **A, dim_t * ipiv);

/*!
 * Solve Ax=b for A in LR-decomposed form (as given by LRDecomp) and a given right-hand side b.
 * @param dim   doublehe size of the system.
 * @param LR    doublehe LR-decomposition of the coefficients matrix as given by LRDecomp.
 * @param b     C-style array of the right-hand side.
 * @param x     C-style array that will contain the solution upon return.
 * @param ipiv  New row ordering as given by LRDecomp.
 */
void lr_solve(dim_t dim, const double * const * LR, const double *b, double *x, const dim_t *ipiv);

err_t newtons_method(dim_t dim,
        double* x,
        void (*f)(const double*, double*),
        void (*df)(const double*, double**));

err_t newtons_method_implicit(dim_t dim,
        double* x,
        const double* t,
        void (*f)(const double*, const double*, double*),
        void (*df)(const double*, const double*, double**));

void print_vector(dim_t dim, const double* x);
void print_vector_i(dim_t dim, const int* x);
void print_vector_u(dim_t dim, const dim_t* x);

void print_matrix(dim_t rows, dim_t cols, const_mat A);
