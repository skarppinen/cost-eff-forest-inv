#ifdef _WIN32
    // This seems to be required in Windows so that PoV is actually put to DLL.
    #define EXPORT_FUNCTION __declspec(dllexport)
#else
    #define EXPORT_FUNCTION
#endif
#include "pov.h"
#include <math.h>
#include <assert.h>
#include "random-sweep.h"

StandData * StandData_alloc(unsigned int nS, unsigned int nT,
        unsigned int nA, unsigned int nI) {
    assert(nS > 0U);
    assert(nT > 0U);
    assert(nA > 0U);
    assert(nI > 0U);

    StandData * stda = malloc(sizeof(StandData));
    double * muprior = calloc(nS * nA, sizeof(double));
    double * sigma2prior = calloc(nS * nA, sizeof(double));
    double * sigma2meas = calloc(nI * nS * nA, sizeof(double));
    double * demands = calloc(nT * nA, sizeof(double));

    stda->nS = nS;
    stda->nT = nT;
    stda->nA = nA;
    stda->nI = nI;
    stda->muprior = muprior;
    stda->sigma2prior = sigma2prior;
    stda->sigma2meas = sigma2meas;
    stda->demands = demands;
    return stda;
}

void StandData_free(StandData * stda) {
    free(stda->muprior);
    free(stda->sigma2prior);
    free(stda->sigma2meas);
    free(stda->demands);
    free(stda);
}

VolumePosterior * VolumePosterior_alloc(unsigned int nS, unsigned int nA) {
    assert(nS > 0U);
    assert(nA > 0U);

    VolumePosterior * vp = malloc(sizeof(VolumePosterior));
    double * muplus = calloc(nS * nA, sizeof(double));
    double * sigma2plus = calloc(nS * nA, sizeof(double));
    double * y = calloc(nS * nA, sizeof(double));

    vp->nS = nS;
    vp->nA = nA;
    vp->muplus = muplus;
    vp->sigma2plus = sigma2plus;
    vp->y = y;
    return vp;
}

void VolumePosterior_free(VolumePosterior * vp) {
    free(vp->muplus);
    free(vp->sigma2plus);
    free(vp->y);
    free(vp);
}

SeedVector * SeedVector_alloc(unsigned int n) {
    assert(n > 0);
    SeedVector * sv = malloc(sizeof(SeedVector));
    unsigned long * arr = calloc(n, sizeof(unsigned long)); 

    sv->arr = arr;
    sv->n = n;
    return sv;
}

void SeedVector_free(SeedVector * sv) {
    free(sv->arr);
    free(sv);
}

static void set_to_zero(double * x, const unsigned int n) {
    for (unsigned int i = 0; i < n; i++)
        x[i] = 0.0;
}

static void u_to_l_tri(double * cm_mat, const unsigned int dim) {
	unsigned int step = 0;
	unsigned int colstart = 0;
	for (unsigned int j = 0; j < dim - 1; ++j) {
		step = dim - 1;
		colstart = j * dim;
		for (unsigned int i = j + 1; i < dim; ++i) {
			cm_mat[colstart + i] = cm_mat[colstart + i + step];
			step += dim - 1;
		}
	}
}

void build_problem(RandomSweepStorage *restrict rs, 
                   const VolumePosterior *restrict vp,
                   const double *restrict demands) {
    double r = 0.0;
    for (unsigned int i = 0; i < vp->nA * rs->nT; i++) {
        r += demands[i] * demands[i]; 
    }
    rs->r = r;

    // Computes `rs.C = -2.0 * muplus * transpose(demands)`.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                vp->nS, rs->nT, vp->nA, 
                -2.0, vp->muplus, vp->nS, 
                demands, rs->nT,
                0.0, rs->C, vp->nS);

    // Computes `rs.Qb`. 
    set_to_zero(rs->Qb, vp->nS * vp->nS);
    for (unsigned int a = 0; a < vp->nA; a++) {
        // Here using striding by `nS + 1` to update diagonal elements of Qb.
        // This adds 2 * Sigma^a to Qb. (Sigma^a is diagonal)
        cblas_daxpy(vp->nS, 
                    2.0, vp->sigma2plus + a * vp->nS, 1,
                    rs->Qb, vp->nS + 1); 
    }
    // This adds sum_{a} 2 * muplus^a * transpose(muplus^a) to Qb.
    // syrk only writes to the upper triangle of Qb here.
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                vp->nS, vp->nA, 
                2.0, vp->muplus, vp->nS,
                1.0, rs->Qb, vp->nS);
    u_to_l_tri(rs->Qb, vp->nS); // Fill lower triangle with what is on upper triangle.
}

void print_StandData(const StandData * stda) {
    int stands[] = {0, 9, 15, 24, 32};
    int assortments[] = {0, 0, 1, 1, 2};
    int inventory_methods[] = {2, 1, 0, 2, 0};

    for (unsigned int i = 0; i < 5; i++) {
        int s = stands[i];
        int a = assortments[i];
        int xi = inventory_methods[i];
        printf("Stand = %d, assortment = %d, inventory = %d:\n",
               s, a, xi); 
        printf("Prior mean = %lf, prior var = %lf, meas var = %lf\n",
                stda->muprior[CMERC(s, a, stda->nS)], 
                stda->sigma2prior[CMERC(s, a, stda->nS)],
                stda->sigma2meas[CMERCS(xi, s, a, stda->nI, stda->nS)]);
    }
}

void lognormal_model_random_posterior(
        VolumePosterior *restrict vp,
        const unsigned int *restrict xI,
        const StandData *restrict stda,
        gsl_rng *restrict rng);

void normal_model_random_posterior(
        VolumePosterior * restrict vp,
        const unsigned int *restrict xI,
        const StandData *restrict stda,
        gsl_rng *restrict rng);

const RandomPosteriorFun RANDOM_POSTERIOR_FUNS[MODEL_CHOICE_COUNT] = {
    [MODEL_CHOICE_NORMAL] = normal_model_random_posterior,
    [MODEL_CHOICE_LOGNORMAL] = lognormal_model_random_posterior,
};

void random_posterior_w_seed(
        VolumePosterior * restrict vp,
        const unsigned int *restrict xI,
        const StandData *restrict stda,
        const ModelChoice model,
        const unsigned long seed) {
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus); 
    gsl_rng_set(rng, seed);
    const RandomPosteriorFun random_posterior = RANDOM_POSTERIOR_FUNS[model];
    random_posterior(vp, xI, stda, rng);
    gsl_rng_free(rng);
}

EXPORT_FUNCTION double PoV(const unsigned int *restrict xI,
           RandomSweepStorage *restrict rs,
           VolumePosterior *restrict vp,
           const StandData *restrict stda, 
           const SeedVector *restrict seeds,  
           const unsigned int maxsweeps,
           const unsigned int ninits,
           const ModelChoice model) {
    if (model < 0 || model >= MODEL_CHOICE_COUNT) {
        printf("Invalid model choice.\n");
        printf("Model choice should be an integer between [%d, %d] (inclusive).\n",
               0, MODEL_CHOICE_COUNT - 1);
        exit(1);
    }
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus); 
    const RandomPosteriorFun random_posterior = RANDOM_POSTERIOR_FUNS[model];
    double out = 0.0, result = 0.0;
    for (unsigned int n = 1; n <= seeds->n; n++) {
        gsl_rng_set(rng, seeds->arr[n - 1]);
        random_posterior(vp, xI, stda, rng);
        build_problem(rs, vp, stda->demands); 
        random_sweep(&result, rs, maxsweeps, ninits, rng); 
        out += (result - out) / n; // Running mean. 
    }
    gsl_rng_free(rng);
    return -out;
}


