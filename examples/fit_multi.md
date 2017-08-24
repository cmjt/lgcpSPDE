### Function fit.multi()

This function was written specifically for the analysis of the British
Trust for Ornithology's (BTO) Garden Bird Feeding Survey (GBFS) data
carried out by Jones-Todd et al. (In Press). The model is given by,

Each **x**<sub>**j**</sub>{**s**<sub>*i*</sub>, *t*}(*j* = 1, 2, 3) is a
spatio-temporal random effect modelled by a SPDE model Lindgren, Rue,
and Lindström (2011), which follows an AR(1) process over time with
parameter *ρ*<sub>*i*</sub>. Each *α*<sub>⋅</sub> is an intercept term
for each component of the model referring to the house sparrows. The
parameters *β*<sub>⋅</sub> and *γ*<sub>⋅</sub> are scaling Blangiardo et
al. (2013 ~Chapter 8) or 'interaction' parameters to the spatio-temporal
random fields of which they are coefficients. That is, each shared
random field (i.e., a random field that appears in more than one linear
predictor) represents the shared inter- or intra-species spatial
auto-correlation over time. Each *β*<sub>⋅</sub> or *γ*<sub>⋅</sub>
parameter represents the magnitude and direction of this spatial
similarity.

This function returns an object (Rue, Martino, and Chopin (2009)) and
takes the following arguments:

-   **locs** A matrix of BTO site locations
-   **mesh** Delauney triangulation of the UK (see Lindgren, Rue, and
    Lindström (2011))
-   **temp** years of GBFS
-   **binary.response** list of length three. each referring to presence
    = 1, absence = 0 of the bird species sparrowhawk, collared dove,
    house sparrow respectively
-   **density.response** list of length three. each referring to
    non-zero density of the bird species sparrowhawk, collared dove,
    house sparrow respectively.
-   **family **a character vector of length two specifying the assumed
    likelihood of each species' response, by default is
    rep(c("binomial","gamma"),3).
-   **control.time** (optional) supplied if the **temp** argument is
    given to fit a spatio-temporal model. This argument controls the
    model and prior put on the hyperparameters of the model for the
    temporal component of the spatio-temporal model. By default this is
    **list(model = 'ar1', param = list(theta = list(prior='pccor1',
    param = c(0, 0.9))))** which is a pc.prior put on the rho
    coefficient of a AR(1) model with P(rho&gt;0)=0.9.
-   **control.inla** a list which controls the fitting procedures INLA
    uses by default this is **list(strategy='gaussian',int.strategy =
    'eb')** for quick and dirty fitting.
-   **hyper** a list (of length 2) of lists of priors for each copy
    parameter. The first list has length 3 specifying the priors on the
    intra species interaction parameters, c(beta\_1, beta\_2, beta\_3)
    (i.e., sparrowhawk, collared dove, house sparrow resp.). the second
    element is of length 4 specifying the priors on the inter species
    interaction parameters in order these refer to the parameters
    beta\_z3, gamma\_z3, beta\_y3, and gamma\_y3 (see model definition).
    By default each is a N(0,10) (i.e.,
    **list(theta=list(prior='normal', param=c(0,10)))**)
-   **control.compute** a list of fit statistics the user wants INLA to
    return. By default this is **list(dic = TRUE, waic = TRUE,cpo =
    TRUE, config = TRUE)**.
-   **spde.new.params** by default this is NULL. If supplied must be a
    named list with components:  - typical standard deviation to use pc
    priors for hyperparams of spde model, **Psig** - prob for sigma of
    pc prior, **rho0** - typical range to use pc priors for hyperparams
    of spde model, and **Prho** - prob for rho of pc prior (see Martins
    et al. (2014))
-   **verbose** Logical if **TRUE** model fit is output to screen.
-   **...** add inla options to speed up computation (i.e., by giving
    starting values from a previous model)

<!-- -->

    library(lgcpSPDE)
    args(fit.multi)

    ## function (locs = NULL, mesh = NULL, temp = NULL, binary.response = NULL, 
    ##     density.response = NULL, family = rep(c("binomial", "gamma"), 
    ##         3), hyper = list(intra = list(beta_1 = list(theta = list(prior = "normal", 
    ##         param = c(0, 10))), beta_2 = list(theta = list(prior = "normal", 
    ##         param = c(0, 10))), beta_3 = list(theta = list(prior = "normal", 
    ##         param = c(0, 10)))), inter = list(beta_z3 = list(theta = list(prior = "normal", 
    ##         param = c(0, 10))), gamma_z3 = list(theta = list(prior = "normal", 
    ##         param = c(0, 10))), beta_y3 = list(theta = list(prior = "normal", 
    ##         param = c(0, 10))), gamma_y3 = list(theta = list(prior = "normal", 
    ##         param = c(0, 10))))), control.time = list(model = "ar1", 
    ##         param = list(theta = list(prior = "pccor1", param = c(0, 
    ##             0.9)))), control.inla = list(strategy = "gaussian", 
    ##         int.strategy = "eb"), control.compute = list(dic = TRUE, 
    ##         waic = TRUE, cpo = TRUE, config = TRUE), spde.new.params = NULL, 
    ##     verbose = FALSE, link = NULL, ...) 
    ## NULL

To fit a basic model using default priors etc, the following code is all
that is required (having specified **locs**, **mesh**, **temp**,
**binary.response**, and **density.response** appropriately).

    fit <- fit.multi(locs = locs, mesh = mesh, temp = temp, ## site locations, mesh, and time indecies
                     binary.response = z.response, ## list of presence/absencse for species
                     density.response = y.response) ## list of 'density' for species

References
==========

Blangiardo, Marta, Michela Cameletti, Gianluca Baio, and Håvard Rue.
2013. “Spatial and Spatio-Temporal Models with R-INLA.” *Spatial and
Spatio-Temporal Epidemiology* 7. Elsevier: 39–55.

Jones-Todd, C. M, B Swallow, J Illian, and M Toms. In Press. “A
Spatio-Temporal Multi-Species Model of a Semi-Continuous Response.”
*Journal of the Royal Statistical Society Series C*.

Lindgren, Finn, Håvard Rue, and Johan Lindström. 2011. “An Explicit Link
Between Gaussian Fields and Gaussian Markov Random Fields: The
Stochastic Partial Differential Equation Approach.” *Journal of the
Royal Statistical Society: Series B (Statistical Methodology)* 73 (4).
Wiley Online Library: 423–98.

Martins, Thiago G, Daniel P Simpson, Andrea Riebler, Håvard Rue, and
Sigrunn H Sørbye. 2014. “Penalising Model Component Complexity: A
Principled, Practical Approach to Constructing Priors.” *arXiv Preprint
arXiv:1403.4630*.

Rue, Håvard, Sara Martino, and Nicolas Chopin. 2009. “Approximate
Bayesian Inference for Latent Gaussian Models by Using Integrated Nested
Laplace Approximations.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 71 (2). Wiley Online Library:
319–92.
