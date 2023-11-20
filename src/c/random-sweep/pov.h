#include "random-sweep.h"

// Get index of (i, j, k) from array A with slices of size `m` * `n`,
// that is, m rows and n columns.
// Slices are the last dimension.
// The naming comes from Column Major Element: Row, Column, Slice
// For example CMERCS(1, 4, 1, 2, 6) returns row index 1 from column index 4 in slice index 1, with each
// slice having 2 rows and 6 columns.
#define CMERCS(i, j, k, m, n) ((k) * (m) * (n) + (j) * (m) + (i))

typedef struct {
    unsigned int nS; // Number of stands.
    unsigned int nT; // Number of time points.
    unsigned int nA; // Number of assortments.
    unsigned int nI; // Number of different inventory methods.
    double * muprior; // Prior volumes, `nS` x `nA` array.
    double * sigma2prior; // Prior volume variances, `nS` x `nA` array.
    double * sigma2meas; // Volume measurement variances, `nI` x `nS` x `nA` array.
    double * demands; // Demands per time and assortment, `nT` x `nA` array.
} StandData;

StandData * StandData_alloc(unsigned int nS, unsigned int nT,
        unsigned int nA, unsigned int nI);
void StandData_free(StandData * stda); 

typedef struct {
    unsigned int nS; // Number of stands.
    unsigned int nA; // Number of assortments.
    double * muplus; // Posterior mean volumes, `nS` x `nA` matrix.
    double * sigma2plus; // Posterior volume variances, `nS` x `nA` matrix.
    double * y; // Temporary for simulated observations, `nS` x `nA` matrix.
} VolumePosterior;

VolumePosterior * VolumePosterior_alloc(unsigned int nS, unsigned int nA);
void VolumePosterior_free(VolumePosterior * vp);

typedef struct {
    unsigned long * arr; // A vector of length `n` with random seeds.
    unsigned int n; // Number of random seeds.
} SeedVector;
SeedVector * SeedVector_alloc(unsigned int n);
void SeedVector_free(SeedVector * sv);

typedef void (*RandomPosteriorFun)(VolumePosterior *, 
            const unsigned int *, const StandData *,
                                   gsl_rng *);

typedef enum {
    MODEL_CHOICE_NORMAL = 0, 
    MODEL_CHOICE_LOGNORMAL,
    MODEL_CHOICE_COUNT
} ModelChoice;

void build_problem(RandomSweepStorage *restrict rs, 
                   const VolumePosterior *restrict vp,
                   const double *restrict demands);

double PoV(const unsigned int *restrict xI,
           RandomSweepStorage *restrict rs,
           VolumePosterior *restrict vp,
           const StandData *restrict stda, 
           const SeedVector *restrict seeds,  
           const unsigned int maxsweeps,
           const unsigned int ninits,
           const ModelChoice);
