#ifndef H_RANDOM_SWEEP
#define H_RANDOM_SWEEP

#ifdef _WIN32
    // This seems to be required in Windows to avoid linking errors with GSL.
    #define WIN32
    #define GSL_DLL
#endif
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h> // This includes cblas.

// Get index (i, j) in a column major array A with m rows.
// The naming comes from Column Major Element: Row, Column
#define CMERC(i, j, m) ((j) * (m) + (i))

enum {
    RANDOM_SWEEP_CONVERGED = 0,
    RANDOM_SWEEP_NOT_CONVERGED = 1
};

typedef enum {
    SWEEP_STRATEGY_SMALLEST_NEGATIVE = 0, 
    SWEEP_STRATEGY_GREATEST_NEGATIVE = 1
} SweepStrategy;

typedef double (*IncrementFinder)(int *, const double *, const unsigned int);

typedef struct {
    double * X; // `nS` x `nT` decision matrix storing current solution.
    double * Xopt; // `nS` x `nT` decision matrix storing best found solution.
    double * C; // `nS` x `nT` matrix containing `c_t` on each column (see below).
    double * Qb; // `nS` x `nS` symmetric matrix in problem formulation (see below). 
    double * QbX; // `nS` x `nT` temporary for storing matrix product `Qb` * `X`.
    double * Qbx; // `nS` x 1 temporary for storing matrix-vector product `Qb` * `x`.
    unsigned int * order; // `nS` x 1 temporary for storing processing order of stands.
    unsigned int * status_lu; // Lookup vector of length `nS` storing current status of each stand. 
    // Temporary for storing increments corresponding to changes to stand harvesting, length `nT + 1`. 
    double * increments; 
    double r; // Constant in objective function.
    unsigned int nS; // Number of stands.
    unsigned int nT; // Number of time points.
    SweepStrategy strategy;
} RandomSweepStorage;

RandomSweepStorage * RandomSweepStorage_alloc(const unsigned int nS,
                                              const unsigned int nT,
                                              const double r,
                                              const SweepStrategy strategy); 


void RandomSweepStorage_free(RandomSweepStorage * rs);

/*
 * Find a feasible solution using a maximum of `max_sweeps` sweeps and `inits` random initialisations
 * to the binary programming problem 
 *      \sum_{t = 1} [ 0.5 * x_t' * Qb * x_t + c_t' x_t ] + r, 
 * where x_t and c_t are the `t`th column of `rs.X` and `rs.C` respectively, and `rs.X` is constrained such
 * that each row contains at most one 1.
 *
 * The return value is an error code, zero means no error, other values correspond to errors.
 * The first argument `objvalue` will point to the best objective function value found.
*/
int random_sweep(double *restrict objvalue, 
                 RandomSweepStorage *restrict rs, 
                 const unsigned int max_sweeps,
                 const unsigned int inits,
                 gsl_rng * rng);

int random_sweep_w_seed(double * out, 
                        RandomSweepStorage * rs,
                        const unsigned int max_sweeps,
                        const unsigned int inits,
                        const unsigned long int seed);

/* 
 * Evaluate \sum_{t = 1} [ 0.5 * x_t' * Qb * x_t + c_t' x_t ] + r, 
 * where x_t and c_t are the `t`th column of `rs.X`and `rs.C` respectively.
*/
double objective_value(const RandomSweepStorage * rs);

/*
 * Generate a random feasible solution to `X` using RNG `rng`.
*/
void random_feasible_solution(double * X, 
                              unsigned int * status_lu,
                              const unsigned int nS,
                              const unsigned int nT,
                              gsl_rng * rng);

/*
 * Evaluate `QbX` = `Qb` * `X`, where `Qb` is an `m` times `m` symmetric matrix,
 * and `X` is a `m` times `n` matrix. The output is an `m` times `n` matrix.
 * All matrices should be in column-major form.
*/
void init_QbX(double *restrict QbX, 
              const double *restrict Qb, 
              const double *restrict X, 
              const unsigned int m, 
              const unsigned int n); 

/*
 * Copy first `n` elements from `src` to `dest`.
*/
void vec_copy(double *restrict dest, const double *restrict src, const unsigned int n);

// Fill `rs.increments` with all possible increments corresponding to all possible status changes, 
// when stand `stand` is currently in status `status`.
// An "increment" is the change in objective function value, when a single stand changes status.
//void fill_increments(RandomSweepStorage * rs, 
//                     const unsigned int stand, 
//                     const unsigned int status); 
void fill_increments(double *restrict increments,
                     const double *restrict Qb,
                     const double *restrict QbX,
                     const double *restrict C,
                     const unsigned int nS,
                     const unsigned int nT,
                     const unsigned int stand,
                     const unsigned int status);

double find_increment(unsigned int *restrict incr_index_ptr, 
                      const IncrementFinder incf,
                      const unsigned int status,
                      const double *restrict increments, 
                      const unsigned int nT); 

// Update the solution given that `stand` changed status from `status` to `newstatus`. 
void update_solution(RandomSweepStorage * rs, 
                     const unsigned int stand, 
                     const unsigned int status, 
                     const unsigned int newstatus); 

// Update `rs.QbX` given that `stand` changed status from `status` to `newstatus`.
// This function exists so that the product `rs.Qb` * `rs.X` need not be recomputed,
// since `rs.X` changes very little when updates to stand statuses happen.
void update_QbX(double *restrict QbX, 
                const double *restrict Qb, 
                const unsigned int nS,
                const unsigned int stand,
                const unsigned int status,
                const unsigned int newstatus);

/* 
 * Shuffle the first `n` elements of `x` using RNG `rng`.
*/
void shuffle(unsigned int * x, const unsigned int n, gsl_rng * rng);

// Return uniformly distributed long unsigned integer in the interval [0, `n`].
// The interval is inclusive in both ends.
unsigned long int rand_zero_upto_n(const unsigned int n, gsl_rng * rng); 

inline void vec_axpy(double *restrict y, const double alpha, 
                     const double *restrict x, const unsigned int n) {
    cblas_daxpy(n, alpha, x, 1, y, 1);
    return;
}

#endif
