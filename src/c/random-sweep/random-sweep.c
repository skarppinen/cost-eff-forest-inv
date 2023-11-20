#include "random-sweep.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h> // memcpy
#include <math.h> // INFINITY

#ifdef DEBUG
unsigned int approx_eq(double x, double y) {
    return fabs(x - y) < 1e-6 ? 1U : 0U;     
}
#endif

RandomSweepStorage * RandomSweepStorage_alloc(const unsigned int nS,
                                              const unsigned int nT,
                                              const double r,
                                              const SweepStrategy strategy) {
    RandomSweepStorage * rs = malloc(sizeof(RandomSweepStorage)); 
    rs->nS = nS;
    rs->nT = nT;
    rs->r = r;
    rs->X = calloc(nS * nT, sizeof(double));
    rs->Xopt = calloc(nS * nT, sizeof(double));
    rs->C = calloc(nS * nT, sizeof(double));
    rs->Qb = calloc(nS * nS, sizeof(double));
    rs->QbX = calloc(nS * nT, sizeof(double));
    rs->Qbx = calloc(nS, sizeof(double));
    rs->order = calloc(nS, sizeof(unsigned int));
    for (unsigned int i = 0; i < nS; i++)
        rs->order[i] = i;
    rs->status_lu = calloc(nS, sizeof(unsigned int));     
    rs->increments = calloc(nT + 1, sizeof(double));
    rs->strategy = strategy;
    return rs;
}

void RandomSweepStorage_free(RandomSweepStorage * rs) {
    free(rs->X);
    free(rs->Xopt);
    free(rs->C);
    free(rs->Qb);
    free(rs->QbX);
    free(rs->Qbx);
    free(rs->order);
    free(rs->status_lu);
    free(rs->increments);
    free(rs);
}

inline void vec_copy(double *restrict dest, const double *restrict src, const unsigned int n) {
    memcpy(dest, src, sizeof(double) * n);
}

void fill_increments(double *restrict increments,
                     const double *restrict Qb,
                     const double *restrict QbX,
                     const double *restrict C,
                     const unsigned int nS,
                     const unsigned int nT,
                     const unsigned int stand,
                     const unsigned int status) {

    const double qss = Qb[CMERC(stand, stand, nS)];
    double status_const = 0.0;
    if (status > 0U) {
        status_const = -C[CMERC(stand, status - 1, nS)] - QbX[CMERC(stand, status - 1, nS)] + 0.5 * qss;
    }
    increments[0] = status_const;
    for (unsigned int to = 1; to < nT + 1; to++) {
        increments[to] = C[CMERC(stand, to - 1, nS)] + 
                         QbX[CMERC(stand, to - 1, nS)] + 
                         0.5 * qss + status_const; 
    }
    increments[status] = 0.0;
}

void update_solution(RandomSweepStorage * rs, 
                     const unsigned int stand, 
                     const unsigned int status, 
                     const unsigned int newstatus) {
    const unsigned int nS = rs->nS;
    if (status == newstatus) {
        return;
    }
    if (status == 0) {
        // `newstatus` is nonzero here.
        rs->X[CMERC(stand, newstatus - 1, nS)] = 1.0;
        rs->status_lu[stand] = newstatus; 
        return;
    }
    if (newstatus == 0) {
        // `status` is nonzero here.
        rs->X[CMERC(stand, status - 1, nS)] = 0.0;
        rs->status_lu[stand] = newstatus; 
        return;
    }
    // Both are nonzero here.
    rs->X[CMERC(stand, newstatus - 1, nS)] = 1.0;
    rs->X[CMERC(stand, status - 1, nS)] = 0.0;
    rs->status_lu[stand] = newstatus; 
    return;
}


inline void init_QbX(double * QbX, 
              const double *restrict Qb,
              const double *restrict X,
              const unsigned int nS,
              const unsigned int nT) {
    cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper,
                nS, nT, 1.0,
                Qb, nS,
                X, nS,
                0.0, QbX, nS);
}

inline double vec_dot(const double *restrict x,
                      const double *restrict y,
                      const unsigned int n) {
   return cblas_ddot(n, x, 1, y, 1); 
}

inline double symm_quadratic_form(const double *restrict A,
                           const double *restrict x,
                           const unsigned int n,
                           double *restrict tmp) {
    // Compute `A * x` to `tmp`.
    cblas_dsymv(CblasColMajor, CblasUpper,
                n, 1.0, A, n,
                x, 1, 
                0.0, tmp, 1);
    return cblas_ddot(n, x, 1, tmp, 1); 
}

double objective_value(const RandomSweepStorage * rs) {
    const unsigned int nS = rs->nS;
    const unsigned int nT = rs->nT;    
    double obj = 0.0;
    for (unsigned int t = 0; t < nT; ++t) {
        obj += 0.5 * symm_quadratic_form(rs->Qb, rs->X + nS * t, nS, rs->Qbx); // x' * Qb * x
        obj += vec_dot(rs->X + nS * t, rs->C + nS * t, nS); // x_t' * c_t
    }
    obj += rs->r;
    return obj;
}

void update_QbX(double *restrict QbX,
                const double *restrict Qb, 
                const unsigned int nS,
                const unsigned int stand, 
                const unsigned int status, 
                const unsigned int newstatus) {
    if (status == newstatus) {
        return;
    }
    const double * Qb_col = Qb + stand * nS; // `stand`th column of `Qb`.
    if (status == 0) {
        // `newstatus` is nonzero here.
        // Add `Qb_col` to column indexed by `newstatus - 1`th in Qb.
        vec_axpy(QbX + (newstatus - 1) * nS, 1.0, Qb_col, nS);
        return;
    }
    if (newstatus == 0) {
        // `status` is nonzero here.
        vec_axpy(QbX + (status - 1) * nS, -1.0, Qb_col, nS);
        return;
    }
    // Both nonzero here.
    vec_axpy(QbX + (newstatus - 1) * nS, 1.0, Qb_col, nS);
    vec_axpy(QbX + (status - 1) * nS, -1.0, Qb_col, nS);
    return;
}


double smallest_negative(int * iptr, const double * x, const unsigned int n) {
    double minimum = INFINITY;
    for (unsigned int i = 0; i < n; i++) {
        if (x[i] < minimum) {
            *iptr = (int) i;
            minimum = x[i];
        }
    }
    return minimum;
}

double greatest_negative(int * iptr, const double * x, const unsigned int n) {
    double greatest = -INFINITY;
    *iptr = -1;
    for (unsigned int i = 0; i < n; i++) {
        if (x[i] < 0.0 && x[i] > greatest) {
            greatest = x[i];
            *iptr = (int) i;
        }
    }
    return greatest; 
}

double find_increment(unsigned int * incr_index_ptr,
                      const IncrementFinder incf,
                      const unsigned int status,
                      const double * increments,
                      const unsigned int nT) {
    int i = 0;
    int * iptr = &i;
    double result = incf(iptr, increments, nT + 1); 
    if (*iptr < 0) {
        result = 0.0;    
        *incr_index_ptr = status;
        return result;
    }
    *incr_index_ptr = (unsigned int) *iptr; 
    return result;
}

unsigned long int rand_zero_upto_n(const unsigned int n, gsl_rng * rng) {
    return gsl_rng_uniform_int(rng, n + 1);
}

void shuffle(unsigned int * x, const unsigned int n, gsl_rng * rng) {
    unsigned long int j = 0;
    unsigned int tmp = 0;
    if (n > 1) {
        for (unsigned int i = n - 1; i >= 1; i--) {
            j = rand_zero_upto_n(i, rng); 
            tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
    }
}

static void vec_set_to_zero(double * x, const unsigned int n) {
    for (unsigned int i = 0; i < n; i++) 
        x[i] = 0.0;
    return;
}

void random_feasible_solution(double * X, 
                              unsigned int * status_lu,
                              const unsigned int nS,
                              const unsigned int nT,
                              gsl_rng * rng) {
    unsigned long int status = 0;
    vec_set_to_zero(X, nS * nT); 
    for (unsigned int i = 0; i < nS; i++) {
        status = rand_zero_upto_n(nT, rng); 
        status_lu[i] = status;
        if (status == 0) {
            continue;
        }
        X[CMERC(i, status - 1, nS)] = 1.0;
    }
    return;
}

#ifdef DEBUG
double rowsum(const double * X, const unsigned int row, 
                     const unsigned int m, const unsigned int n) {
    double sum = 0.0;
    for (unsigned int j = 0; j < n; j++)
        sum += X[CMERC(row, j, m)]; 
    return sum;
}

unsigned int isfeasible(const double * X, const unsigned int nS,
                        const unsigned int nT) {
    double s = 0.0;
    for (unsigned int i = 0; i < nS; i++) {
        s = rowsum(X, i, nS, nT); 
        if (s < 0.0 || s > 1.0) {
            return 0U;
        }
    }
    return 1U;
}
#endif 

int random_sweep_w_seed(double * out, 
                        RandomSweepStorage * rs,
                        const unsigned int max_sweeps,
                        const unsigned int inits,
                        const unsigned long int seed) {
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    int res = random_sweep(out, rs, max_sweeps, inits, rng); 
    gsl_rng_free(rng);
    return res;
}

const IncrementFinder INCREMENT_FINDERS[2] = {
    [SWEEP_STRATEGY_SMALLEST_NEGATIVE] = &smallest_negative,
    [SWEEP_STRATEGY_GREATEST_NEGATIVE] = &greatest_negative
};

int random_sweep(double * out, 
                 RandomSweepStorage * rs, 
                 const unsigned int max_sweeps, 
                 const unsigned int inits,
                 gsl_rng * rng) {

    assert(inits > 0U);
    assert(max_sweeps > 0U);
    const IncrementFinder incf = INCREMENT_FINDERS[rs->strategy];
    const unsigned int nS = rs->nS;
    const unsigned int nT = rs->nT;
    const unsigned int N = nS * nT;
    unsigned int changed = 0U;
    unsigned int incr_index = 0U; 
    unsigned int stand = 0U;
    unsigned int newstatus = 0U;
    unsigned int status = 0U;
    double incr = 0.0;
    double objvalue = 0.0; 
    double best_objvalue = INFINITY;
    *out = best_objvalue;
    int ret = RANDOM_SWEEP_NOT_CONVERGED;
    // Required to ensure reproducible results wrt RNG.
    for (unsigned int i = 0U; i < nS; i++)
        rs->order[i] = i;
    
    for (unsigned int k = 0; k < inits; k++) {
        random_feasible_solution(rs->X, rs->status_lu, nS, nT, rng); 
        init_QbX(rs->QbX, rs->Qb, rs->X, nS, nT); 
        objvalue = objective_value(rs); 
        if (objvalue < best_objvalue) {
           best_objvalue = objvalue; 
           vec_copy(rs->Xopt, rs->X, N); 
        }
        for (unsigned int i = 0U; i < max_sweeps; i++) {
            shuffle(rs->order, nS, rng); 
            changed = 0U;
            for (unsigned int j = 0U; j < nS; j++) {
                stand = rs->order[j];
                status = rs->status_lu[stand];
                fill_increments(rs->increments, rs->Qb, rs->QbX, rs->C, nS, nT, stand, status); 
                incr = find_increment(&incr_index, incf, status, rs->increments, nT); 
                newstatus = incr_index;
                if (newstatus != status) {
                    changed = 1U;
                    update_solution(rs, stand, status, newstatus); 
                    update_QbX(rs->QbX, rs->Qb, nS, stand, status, newstatus);
                    objvalue += incr; 
#ifdef DEBUGVERBOSE
                    if (incr < 0.0) {
                        printf("init = %u, sweep = %u, stand = %u, objective updates to %lf\n", 
                                k, i, j, objvalue);
                    }
#endif
                    if (objvalue < best_objvalue) {
                        best_objvalue = objvalue;  
                        vec_copy(rs->Xopt, rs->X, N); 
                    }
                }
#ifdef DEBUG
                assert(incr <= 0.0);
                assert(isfeasible(rs->X, rs->nS, rs->nT));
                assert(approx_eq(objective_value(rs), objvalue)); 
#endif
            }
            if (!changed) { // Early exit from sweeping loop.
#ifdef DEBUGVERBOSE
                printf("EARLY EXIT\n");
#endif
                ret = RANDOM_SWEEP_CONVERGED;
                break;
            }
        }
    }
    *out = best_objvalue;
    return ret;
}
