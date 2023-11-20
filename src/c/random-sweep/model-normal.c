#include "pov.h"
#include <math.h>

static void normal_posterior_mean_var(double * mupost, double * sigma2post,
                                      double y, double muprior, double sigma2prior,
                                      double sigma2meas) { 
    *sigma2post = 1.0 / ((1.0 / sigma2prior) + (1.0 / sigma2meas));
    *mupost = *sigma2post * (muprior / sigma2prior + y / sigma2meas);
}

static void normal_model_simulate_data(double *restrict y, const unsigned int *restrict xI, 
                   const StandData *restrict stda, gsl_rng *restrict rng) {
    double priorm = 0.0;
    double priorv = 0.0;
    double measv = 0.0;
    for (unsigned int a = 0; a < stda->nA; a++) {
        for (unsigned int s = 0; s < stda->nS; s++) {
            measv = stda->sigma2meas[CMERCS(xI[s], s, a, stda->nI, stda->nS)]; 
            priorm = stda->muprior[CMERC(s, a, stda->nS)];
            priorv = stda->sigma2prior[CMERC(s, a, stda->nS)];
            y[CMERC(s, a, stda->nS)] = priorm + gsl_ran_gaussian_ziggurat(rng, sqrt(priorv + measv)); 
        }
    }
}

void normal_model_compute_posterior(VolumePosterior *restrict vp,
                       const unsigned int *restrict xI,
                       const StandData *restrict stda) {
    double priorm = 0.0;
    double priorv = 0.0;
    double measv = 0.0;
    for (unsigned int a = 0; a < stda->nA; a++) {
        for (unsigned int s = 0; s < stda->nS; s++) {
            measv = stda->sigma2meas[CMERCS(xI[s], s, a, stda->nI, stda->nS)]; 
            priorm = stda->muprior[CMERC(s, a, stda->nS)];
            priorv = stda->sigma2prior[CMERC(s, a, stda->nS)];
            normal_posterior_mean_var(vp->muplus + CMERC(s, a, vp->nS), 
                                      vp->sigma2plus + CMERC(s, a, vp->nS),
                                      vp->y[CMERC(s, a, vp->nS)], priorm, priorv, measv);
        }
    }
}

void normal_model_random_posterior(VolumePosterior *restrict vp,
                                   const unsigned int *restrict xI,
                                   const StandData *restrict stda,
                                   gsl_rng *restrict rng) {
    normal_model_simulate_data(vp->y, xI, stda, rng);
    normal_model_compute_posterior(vp, xI, stda);
}
