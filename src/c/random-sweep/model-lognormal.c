#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include "pov.h"

// Setting for numerical integration.
#define WORKSPACE_SIZE (2000)
#define INTEGRAL_RELTOL (1e-4)
#define INTEGRAL_ABSTOL (0.0)

typedef struct {
    double muprior;
    double sigma2prior;
    double sigma2meas; 
    double y;

    // Save these also to avoid recomputing. 
    double std_mu;
    double std_sigma; 
    double nu_p;
    double sqrt_lam_p;
    double postmean;
    double lognormconst;
} LognormalParams;

static void LognormalParams_print(const LognormalParams * p) {
    printf("LognormalParams { muprior = %lf, sigma2prior = %lf, sigma2meas = %lf, y = %lf } \n"
           "additional fields: std_mu = %lf, std_sigma = %lf, nu_p = %lf, sqrt_lam_p = %lf\n",
            p->muprior, p->sigma2prior, p->sigma2meas, p->y, p->std_mu, p->std_sigma, 
            p->nu_p, p->sqrt_lam_p); 
}

// Unnormalised logposterior density log(p(v \mid y))
static double lp_unnorm_v(double v, const LognormalParams * p) {
    const double muprior = p->muprior;
    //const double sigma2prior = p->sigma2prior;
    const double sigma2meas = p->sigma2meas;
    const double y = p->y;

    const double vsq = v * v;
    const double mupriorsq = muprior * muprior;
    //const double nu_p = 2.0 * log(muprior) - 0.5 * log(mupriorsq + sigma2prior);
    const double nu_m = 2.0 * log(v) - 0.5 * log(vsq + sigma2meas);
    //const double lam_p = log1p(sigma2prior / mupriorsq);
    const double lam_m = log1p(sigma2meas / vsq);
    if (lam_m <= 0.0) return -INFINITY;
    const double logyminusnu_div_sqrt_lam_m = (log(y) - nu_m) / sqrt(lam_m); 
    const double logvminusnup_div_sqrt_lam_p = (log(v) - p->nu_p) / p->sqrt_lam_p;

    const double lout = -0.5 * log(lam_m) - 
        0.5 * logyminusnu_div_sqrt_lam_m * logyminusnu_div_sqrt_lam_m  -
        log(v) - 0.5 * logvminusnup_div_sqrt_lam_p * logvminusnup_div_sqrt_lam_p;

    if (isnan(lout)) {
        printf("Error in numerical integration: lout is NaN!\n");
        LognormalParams_print(p);
        printf("v = %lf\n", v);
        printf("vsq = %lf\n", vsq);
        printf("mupriorsq = %lf\n", mupriorsq);
        printf("nu_p = %lf\n", p->nu_p);
        printf("nu_m = %lf\n", nu_m);
        printf("lam_p = %lf\n", p->sqrt_lam_p * p->sqrt_lam_p);
        printf("lam_m = %lf\n", lam_m);
        printf("Exiting...\n");
        exit(1);
    }
    if (lout > 0.0 && isinf(lout)) {
        printf("Error in numerical integration: lout is positive Inf!\n");
        LognormalParams_print(p);
        printf("v = %lf\n", v);
        printf("vsq = %lf\n", vsq);
        printf("mupriorsq = %lf\n", mupriorsq);
        printf("nu_p = %lf\n", p->nu_p);
        printf("nu_m = %lf\n", nu_m);
        printf("lam_p = %lf\n", p->sqrt_lam_p * p->sqrt_lam_p);
        printf("lam_m = %lf\n", lam_m);
        printf("Exiting...\n");
        exit(1);
    }
    return lout;
}

/*
 * Standardised versions of normalising constant and unnormalised moment integrals.
 */
static double p_unnorm_vstd(double vstd, void * params) {
    const LognormalParams * p = (LognormalParams *) params; 
    const double v = p->std_mu + p->std_sigma * vstd;
    return exp(log(p->std_sigma) + lp_unnorm_v(v, p));
}
static double unnorm_mom1_integrand_vstd(double vstd, void * params) {
    const LognormalParams * p = (LognormalParams *) params;
    const double v = p->std_mu + p->std_sigma * vstd;
    return exp(log(p->std_sigma) + log(v) + lp_unnorm_v(v, p));
}

static double norm_var_integrand_vstd(double vstd, void * params) {
    const LognormalParams * p = (LognormalParams *) params;
    const double v = p->std_mu + p->std_sigma * vstd;
    return exp(log(p->std_sigma) + 2.0 * log(fabs(v - p->postmean)) - p->lognormconst + lp_unnorm_v(v, p));
}

void lognormal_model_simulate_data(double *restrict y, 
        const unsigned int *restrict xI,
        const StandData *restrict stda, 
        gsl_rng *restrict rng) {

    for (unsigned int a = 0; a < stda->nA; a++) {
        for (unsigned int s = 0; s < stda->nS; s++) {
            const double measv = stda->sigma2meas[CMERCS(xI[s], s, a, stda->nI, stda->nS)]; 
            const double priorm = stda->muprior[CMERC(s, a, stda->nS)];
            const double priorv = stda->sigma2prior[CMERC(s, a, stda->nS)];
            const double v = gsl_ran_lognormal(rng, 
                    2.0 * log(priorm) - 0.5 * log(priorm * priorm + priorv),
                    sqrt(log1p(priorv / (priorm * priorm))));

            y[CMERC(s, a, stda->nS)] = gsl_ran_lognormal(rng, 
                    2.0 * log(v) - 0.5 * log(v * v + measv),
                    sqrt(log1p(measv / (v * v))));
        }
    }

}

static void integrate_mean_var(double * mean, double * var, LognormalParams * p, 
                        gsl_integration_workspace * ws) {
    double abserr = 0.0;
    p->std_sigma = sqrt(p->sigma2prior);
    p->std_mu = p->muprior; 
    const double lower = -p->std_mu / p->std_sigma; // Lower integration bound for standardised variable. 
    const double upper = 100.0; // Upper integration bound for standardised variable (practically zero after this).

    // All of the integrals below need these parameters.
    const double mupriorsq = p->muprior * p->muprior;
    p->nu_p = 2.0 * log(p->muprior) - 0.5 * log(mupriorsq + p->sigma2prior);
    p->sqrt_lam_p = sqrt(log1p(p->sigma2prior / mupriorsq));

    // Find normalising constant.
    gsl_function I_norm_const;
    I_norm_const.function = &p_unnorm_vstd;
    I_norm_const.params = p; 
    double normconst = 0.0;
    gsl_integration_qag(&I_norm_const, lower, upper, INTEGRAL_ABSTOL, INTEGRAL_RELTOL, 
            WORKSPACE_SIZE, GSL_INTEG_GAUSS61, ws, &normconst, &abserr); 
    if (isnan(normconst)) {
        printf("Error in numerical integration: normalisation constant is NaN\n");
        printf("Parameters causing trouble were:\n");
        LognormalParams_print(p);
        exit(1);
    }
    if (fabs(normconst) <= 0.0) {
        printf("Error in numerical integration: normalisation constant is zero\n");
        printf("Parameters causing trouble were:\n");
        LognormalParams_print(p);
        exit(1);
    }

    // Find "unnormalised" first moment.
    gsl_function I_unnorm_mom1;
    I_unnorm_mom1.function = &unnorm_mom1_integrand_vstd;
    I_unnorm_mom1.params = p;
    double unnorm_mom1 = 0.0;
    gsl_integration_qag(&I_unnorm_mom1, lower, upper, INTEGRAL_ABSTOL, INTEGRAL_RELTOL, 
            WORKSPACE_SIZE, GSL_INTEG_GAUSS61, ws, &unnorm_mom1, &abserr); 
    if (isnan(unnorm_mom1)) {
        printf("Error in numerical integration: unnormalised first moment is NaN\n"); 
        printf("Parameters causing trouble were:\n");
        LognormalParams_print(p);
        exit(1);
    }
    if (unnorm_mom1 <= 0.0) {
        printf("Error in numerical integration: unnormalised first moment is zero\n");
        printf("Parameters causing trouble were:\n");
        LognormalParams_print(p);
        exit(1);
    }
    const double m = unnorm_mom1 / normconst; 
    if (m <= 0.0) {
        printf("Error in numerical integration: got negative posterior mean!\n");
        printf("Parameters causing trouble were:\n");
        LognormalParams_print(p);
        exit(1);
    }
    p->postmean = m;
    p->lognormconst = log(normconst); 

    gsl_function I_norm_var;
    I_norm_var.function = &norm_var_integrand_vstd;
    I_norm_var.params = p;
    double norm_var = 0.0;
    gsl_integration_qag(&I_norm_var, lower, upper, INTEGRAL_ABSTOL, INTEGRAL_RELTOL, 
            WORKSPACE_SIZE, GSL_INTEG_GAUSS61, ws, &norm_var, &abserr); 
    if (isnan(norm_var)) {
        printf("Error in numerical integraion: posterior variance is NaN\n"); 
        printf("Parameters causing trouble were:\n");
        LognormalParams_print(p);
        exit(1);
    }
    if (norm_var <= 0.0) {
        printf("Error in numerical integration: got negative posterior variance\n"); 
        printf("Parameters causing trouble were:\n");
        LognormalParams_print(p);
        exit(1);
    }
    
    *mean = m; 
    *var = norm_var;
    
}

void lognormal_model_approx_posterior(VolumePosterior * restrict vp,
        const unsigned int *restrict xI,
        const StandData *restrict stda) {
    gsl_integration_workspace * ws = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
    for (unsigned int a = 0; a < stda->nA; a++) {
        for (unsigned int s = 0; s < stda->nS; s++) {
            LognormalParams p = {0};
            p.sigma2meas = stda->sigma2meas[CMERCS(xI[s], s, a, stda->nI, stda->nS)]; 
            p.muprior = stda->muprior[CMERC(s, a, stda->nS)];
            p.sigma2prior = stda->sigma2prior[CMERC(s, a, stda->nS)];
            p.y = vp->y[CMERC(s, a, vp->nS)]; 
            integrate_mean_var(vp->muplus + CMERC(s, a, vp->nS),
                               vp->sigma2plus + CMERC(s, a, vp->nS), &p, ws);
        }
    }
    gsl_integration_workspace_free(ws);
}

void lognormal_model_random_posterior(VolumePosterior * restrict vp,
        const unsigned int *restrict xI,
        const StandData *restrict stda,
        gsl_rng * rng) {
    lognormal_model_simulate_data(vp->y, xI, stda, rng);
    lognormal_model_approx_posterior(vp, xI, stda);
}




