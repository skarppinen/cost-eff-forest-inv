using RCall

function R_init_packages()
    R"""
        suppressPackageStartupMessages(library(dplyr))
        suppressPackageStartupMessages(library(tidyr))
        suppressPackageStartupMessages(library(here))
    """
end

function R_get_stand_data(; id::Integer = 3)
    R"""
    suppressPackageStartupMessages(library(here))
    filename <- paste0("test-data-", $(string(id)), ".rds")
    filepath <- file.path(here(), "no-vc", "data", "clean", filename)
    o <- readRDS(filepath)
    """
    muprior = rcopy(R"""o$muprior""");
    varprior = rcopy(R"""o$sdprior""") .^ 2.0;
    measvars = rcopy(R"""o$meas_sds""") .^ 2.0;
    demands = rcopy(R"""o$demands""");
    costs = rcopy(R"""o$costs""");
    stand_ids = rcopy(R"""o$stand_ids""");
    areas = rcopy(R"""o$areas""");
    coordinates = rcopy(R"""o$coordinates""");
    (; muprior, varprior, measvars, demands, costs, stand_ids, areas, coordinates);
end

function R_get_prior_data_and_demands()
    R_init_packages()
	mu0_df = rcopy(R"""
        fp <- file.path(here(), "no-vc", "data", "clean", "volume-data-2.rds")
		o <- readRDS(fp)
		mu0_df <- pivot_wider(select(o$priors, -stand, -prior_sd), names_from = "assortment",
			values_from = "prior_mean")
		mu0_df
	""");
	sigma0_df = rcopy(R"""
			fp <- file.path(here(), "no-vc", "data", "clean", "volume-data-2.rds")
			o <- readRDS(fp)
			sigma0_df <- pivot_wider(select(o$priors, -stand, -prior_mean), names_from = "assortment", values_from = "prior_sd")
			sigma0_df
	""");
	demands_df = rcopy(R"""
			fp <- file.path(here(), "no-vc", "data", "clean", "volume-data-2.rds")
			o <- readRDS(fp)
			demands_df <- o$demands
			demands_df
    """);
	assortments = [:pine, :spruce, :deciduous];
	mu0s = Matrix(mu0_df[:, assortments]);
	sigma0s = Matrix(sigma0_df[:, assortments]);
	D = Matrix(demands_df[:, assortments]);
	
	mu0s, sigma0s, D;
end

function R_get_meas_sds()
    R_init_packages();
	meas_sds = rcopy(R"""
		fp <- file.path(here(), "no-vc", "data", "clean", "volume-data-2.rds")
		o <- readRDS(fp)
		meas_sd_df <- pivot_wider(select(o$meas_sd, -stand), names_from = "assortment",
			values_from = "meas_sd")
		assortments <- c("pine", "spruce", "deciduous")
		mat <- as.matrix(meas_sd_df[, assortments])
		meas_sds <- array(mat, c(3, 200, 3), dimnames = list(c(1, 2, 3), NULL, assortments)) 
		aperm(meas_sds, c(2, 1, 3)) # Dimensions: Stand, inventory, assortment
		""");
	meas_sds;
end
