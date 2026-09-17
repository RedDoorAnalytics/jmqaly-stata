{smcl}
{* *! version 1.6.5  17sep2026}{...}
{vieweralsosee "merlin" "help merlin"}{...}
{vieweralsosee "merlin postestimation" "help merlin_postestimation"}{...}
{viewerjumpto "Syntax" "jmqaly##syntax"}{...}
{viewerjumpto "Description" "jmqaly##description"}{...}
{viewerjumpto "Options" "jmqaly##options"}{...}
{viewerjumpto "Remarks" "jmqaly##remarks"}{...}
{viewerjumpto "Examples" "jmqaly##examples"}{...}
{viewerjumpto "References" "jmqaly##references"}{...}
{viewerjumpto "Author" "jmqaly##author"}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:jmqaly}
{it:newvarname}
{cmd:,}
{cmd:{c -(}}{opt timev:ar(varname)}{cmd:|}{opt life:time}{cmd:{c )-}}
[{it:options}]

{pstd}
{cmd:jmqaly} is for use after {helpb merlin} has fitted a joint
longitudinal-survival model. It is not an estimation command. It requires
StataNow {bf:19.5} and {cmd:merlin} {bf:2.5.0} specifically.

{synoptset 30 tabbed}{...}
{synopthdr:options}
{synoptline}
{syntab:Main}
{synopt:{opt timev:ar(varname)}}time points at which to evaluate the
prediction; specify this or {opt lifetime}{p_end}
{synopt:{opt life:time}}integrate to the whole horizon rather than a time grid;
specify this or {opt timevar()}{p_end}
{synopt:{opt at1(varname # [...])}}covariate pattern for the first
prediction{p_end}
{synopt:{opt at2(varname # [...])}}covariate pattern for the second
prediction; requesting it also produces the difference{p_end}
{synopt:{opt qolm:odel(#)}}which submodel holds quality of life; required when
the fit has more than two submodels{p_end}
{synopt:{opt ci}}also produce 95% confidence intervals{p_end}

{syntab:What to predict}
{synopt:{opt surv:ival}}predict the survival probability instead of the
QALY{p_end}
{synopt:{opt rmst}}predict the restricted mean survival time instead of the
QALY{p_end}
{synopt:{opt drate(#)}}annual discount rate applied to the QALY; default
{cmd:drate(0)}{p_end}
{synopt:{opt dme:thod(method)}}{cmd:continuous} (the default) or
{cmd:discrete}{p_end}
{synopt:{opt tol:erance(#)}}with {opt lifetime}: the remaining contribution to
stop at, as a fraction of the total; default {cmd:tolerance(1e-4)}{p_end}
{synopt:{opt maxh:orizon(#)}}with {opt lifetime}: hard cap on the horizon
searched; default 500 times the largest observed event time{p_end}

{syntab:Random effects}
{synopt:{opt fixedonly}}predict with the random effects set to zero, rather
than integrating over them{p_end}
{synopt:{opt stand:ardise}}average over the observed covariate distribution
rather than setting other covariates to zero; see {help jmqaly##standardise:below}{p_end}

{syntab:Quality of life scale}
{synopt:{opt transform}}quality of life was modelled on the logit scale;
back-transform it before integrating{p_end}
{synopt:{opt frommin(#)}}lower end of the back-transformed range; default
{cmd:frommin(0)}{p_end}
{synopt:{opt frommax(#)}}upper end of the back-transformed range; default
{cmd:frommax(1)}{p_end}
{synopt:{opt tomin(#)}}lower end of the utility range to rescale onto;
default {cmd:tomin(0)}{p_end}
{synopt:{opt tomax(#)}}upper end of the utility range to rescale onto;
default {cmd:tomax(1)}{p_end}

{syntab:Numerical}
{synopt:{opt qintpoints(#)}}Gauss-Legendre nodes for the integral over time;
default {cmd:qintpoints(30)}{p_end}
{synopt:{opt iter:ate(#)}}passed to {helpb predictnl}; default
{cmd:iterate(100)}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:jmqaly} calculates restricted quality-adjusted life years (QALYs)
following a joint longitudinal-survival model fitted with {helpb merlin}. For
a horizon {it:t} it computes

{p 12 12 2}
QALY({it:t}) = ∫₀ᵗ Q({it:u}) S({it:u}) exp(-{it:r} {it:u}) d{it:u}

{pstd}
where Q({it:u}) is the mean quality of life trajectory, S({it:u}) the survival
probability, and {it:r} the discount rate from {opt drate()}. The integral is
evaluated by Gauss-Legendre quadrature.

{pstd}
By default the prediction is {it:marginal} with respect to the random effects:
Q and S are multiplied at each quadrature node of the random effect
distribution, and {cmd:merlin} integrates the product over that distribution.
The estimand is therefore E[Q({it:u})S({it:u})], which carries Cov(Q,S) -- the
association between quality of life and survival that the joint model
estimates.

{pstd}
{cmd:jmqaly} creates {it:newvarname}{cmd:_at1}, and with {opt at2()} also
{it:newvarname}{cmd:_at2} and {it:newvarname}{cmd:_diff}. With {opt ci} each of
those gains {cmd:_lci} and {cmd:_uci}. Confidence intervals come from the
numerical delta method, via {helpb predictnl}.

{pstd}
The survival function is resolved from the family the survival submodel was
fitted with, so any survival family {cmd:merlin} can predict from may be used.
There must be exactly one survival submodel, and the fit must be two-level.

{pstd}
The quality of life trajectory is taken on its natural scale, with the
submodel's inverse link applied, so {cmd:gaussian}, {cmd:bernoulli},
{cmd:gamma}, {cmd:poisson} and {cmd:beta} submodels all work. {cmd:merlin}
supplies no expected value for its ordinal families, so those are refused by
name.

{pstd}
Cross-submodel association structures are refused: {cmd:EV[]}, {cmd:dEV[]},
{cmd:iEV[]}, {cmd:d2EV[]} and {cmd:XB[]}. {cmd:jmqaly} evaluates the survival
submodel at quadrature nodes of its own choosing, and these make the hazard at
each node depend on another submodel evaluated at the same point, which that
evaluation does not carry. Share the random effects between submodels
instead.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt timevar(varname)} gives the time points at which the prediction is
evaluated, producing a QALY as a function of the horizon. At {it:t}=0 the QALY
and the RMST are zero.

{phang}
{opt lifetime} integrates to the whole horizon instead, producing one number
rather than a curve. Exactly one of {opt timevar()} and {opt lifetime} must be
given: neither is a default, because a QALY curve and a lifetime QALY are
different estimands and the caller should say which is wanted. See
{help jmqaly##lifetime:Remarks on lifetime QALYs} for what "the whole horizon"
resolves to and when it does not exist.

{phang}
{opt at1(varname # [...])} and {opt at2(varname # [...])} fix covariates at
the values given, exactly as {opt at()} does in {helpb merlin_postestimation:predict}
after {cmd:merlin}. Covariates not listed are set to zero unless
{opt standardise} is specified. Requesting {opt at2()} also produces the
difference {cmd:_diff}, formed inside a single {helpb predictnl} expression so
that its confidence interval accounts for the covariance between the two
predictions.

{phang}
{opt qolmodel(#)} names the submodel holding quality of life. With two
submodels it is inferred -- the one that is not the survival submodel -- and
with more than two it is required rather than guessed at.

{dlgtab:What to predict}

{phang}
{opt survival} and {opt rmst} replace the QALY with the survival probability
or the restricted mean survival time for the same covariate pattern. They are
the two components of the QALY integrand, and are useful for checking a fit.
Only one may be given, and neither may be combined with {opt drate()} or
{opt transform}, which apply only to the QALY.

{pstd}
They are diagnostics for the QALY rather than general-purpose predictions, so
a fit {cmd:jmqaly} cannot compute a QALY from is refused whichever of the three
is asked for. For a marginal survival curve from such a fit, use
{helpb merlin_postestimation:predict} after {cmd:merlin}.

{phang}
{opt drate(#)} discounts future quality-adjusted survival, weighting the
integrand by a discount factor.

{phang}
{opt dmethod(method)} selects the factor. {cmd:continuous}, the default,
weights by exp(-{it:r} {it:u}). {cmd:discrete} weights by
(1+{it:r})⁻ᵘ, which is what several health-technology guidelines
specify. The two differ by about 0.6% at ten years at {it:r}=0.035, the
discrete factor discounting slightly less.

{dlgtab:Random effects}

{phang}
{opt fixedonly} sets the random effects to zero instead of integrating over
them. The result is the QALY for a subject at the centre of the random effect
distribution, which is not the population average whenever Cov(Q,S) is
non-zero.

{marker standardise}{...}
{phang}
{opt standardise} averages over the observed distribution of the covariates
not named in {opt at1()} or {opt at2()}, instead of setting them to zero. It
requires
{cmd:merlin} support for standardisation after a model with random effects;
until then {cmd:merlin} declines it, since a joint model always has them.

{dlgtab:Quality of life scale}

{phang}
{opt transform} covers the case where a continuous utility was modelled on the
logit scale -- which is worth doing, since it keeps predicted utilities in a
plausible range. The linear predictor is passed through the inverse logit and
then mapped linearly from [{opt frommin()}, {opt frommax()}] onto
[{opt tomin()}, {opt tomax()}]. The defaults map (0,1) onto (0,1), i.e. the
inverse logit alone; to target the EQ-5D range use
{cmd:transform tomin(-0.594) tomax(1)}.

{phang}
{opt transform} is only needed when the link is the identity. For a submodel
fitted with a non-identity link -- {cmd:family(bernoulli)}, say -- the inverse
link is applied already, and {opt transform} should not be given.

{dlgtab:Numerical}

{phang}
{opt qintpoints(#)} sets the number of Gauss-Legendre nodes for the integral
over time. The default of 30 is converged: against 200 nodes it differs by
about 2e-08 relative, and the relative change between 30 and 100 nodes is of
the same order. Raising it costs time roughly in proportion when {opt ci} is
requested, since every node is re-evaluated for each perturbed parameter.

{phang}
{opt tolerance(#)} is used with {opt lifetime} and sets how much of the
integral may remain unaccounted for when the search stops, as a fraction of
the total. The default is {cmd:tolerance(1e-4)}. It is a bound on the
{it:remaining contribution}, not on survival: under a shared random effect the
marginal survival has a much heavier tail than any individual's, so survival
can sit at 1e-4 over a long stretch while contributing almost nothing to the
QALY. Bounding the contribution measures the quantity actually wanted.

{phang}
{opt maxhorizon(#)} is used with {opt lifetime} and caps how far the search
will go, defaulting to 500 times the largest observed event time. Reaching it
is a refusal rather than a silent truncation.

{marker lifetime}{...}
{title:Remarks on lifetime QALYs}

{pstd}
{opt lifetime} does not compute an integral to infinity. It doubles the horizon
until the accumulated QALY stops moving, and then reports
{it:the horizon it used} along with the contribution it estimates is still
missing. Both are
echoed and returned in {cmd:r(horizon)} and {cmd:r(remaining)}. "Integrated to
232, with an estimated 3e-05 of the total remaining" is a statement a reader
can check; "integrated to infinity" is not.

{pstd}
The search works on the estimand rather than on survival. Writing {it:d_k} for
the increase in the QALY over the {it:k}th doubling, the ratio
{it:rho} = {it:d_k}/{it:d_k-1} measures how fast the integrand is decaying, and
while {it:rho} < 1 the contribution still to come is estimated as
{it:d_k} {it:rho}/(1-{it:rho}). For a survival tail behaving like
{it:t}^-{it:a} this ratio tends to 2^(1-{it:a}), which is below 1 exactly when
{it:a} > 1 -- the condition under which the integral converges. The test is
therefore the convergence condition itself rather than a proxy for it.

{pstd}
{bf:The integral does not always exist.} A model with a surviving fraction --
a Gompertz with negative shape, for instance -- leaves an integrand that never
vanishes, and undiscounted its lifetime QALY is unbounded. {cmd:jmqaly} refuses
rather than returning a number, and says which of the possible causes it
found: increments that are not shrinking (the integral diverges; use
{opt drate()}), or a {opt tolerance()} below what the quadrature can
demonstrate (loosen it, or raise {opt qintpoints()} to lower the noise floor,
which does not remove it).

{pstd}
{bf:A negative utility trajectory is a warning, not a refusal.} If the fitted
quality-of-life trajectory crosses zero within the horizon, the lifetime QALY
accumulates negative utility from that point. With discounting the integral
still converges, so the number is well defined and {cmd:jmqaly} returns it --
but it names the crossing point, because a utility extrapolated far beyond the
data is a modelling question rather than a numerical one. {opt transform}
bounds the utility scale and removes the issue.

{pstd}
{bf:Two assumptions are worth stating.} The remaining-contribution estimate
assumes increments that fall and then stop falling; an integrand whose
increments fall and then grow without bound would be read as having settled.
Bounded utility rules that out, and {opt transform} bounds it. And a lifetime
horizon extrapolates the fitted model far past the data, including its
survival: no background-mortality anchoring is applied, so the estimand is
model-only and not all-cause.

{marker remarks}{...}
{title:Remarks}

{pstd}
The number of quadrature points used for the random effects is
{cmd:merlin}'s, not {cmd:jmqaly}'s.

{pstd}
{opt survival} is 1 at {it:t}=0, with degenerate confidence limits, matching
{cmd:merlin}'s own {helpb merlin_postestimation:predict}. For the families
whose cumulative hazard is an integral -- {cmd:weibull}, {cmd:exponential},
{cmd:gompertz} and the rest -- that is a theorem, the integral being over an
empty interval. For {cmd:family(rp)} and {cmd:family(logchazard)} the
cumulative hazard is exp(spline(log {it:t})) instead, and the value at zero
follows the sign of the spline's lower-tail slope: positive gives S(0)=1,
negative gives S(0)=0. 1 is correct under {cmd:merlin}'s default knot
placement, which is what both packages assume, and under it a negative slope
is not reachable by accident -- though a deliberate {opt knots()} with the
lower boundary far below the data can produce one. No runtime check is
attempted, because the coefficient is recoverable from {cmd:e()} only for a
baseline-only {cmd:rp} model: with a time-dependent term or random effects the
limit is per-observation, so a check that passed would be silently wrong
outside that one case.

{pstd}
The QALY and the RMST are zero at {it:t}=0 for every family, for the same
reason as the first case: an integral over an interval of zero length is zero
whatever the integrand.

{marker ltruncated}{...}
{pstd}
{bf:Delayed entry.} If the survival submodel was fitted with
{opt ltruncated()}, {cmd:jmqaly} still accrues from time zero: the estimand is
the QALY from the origin, not the QALY from each subject's entry time, and the
survival it integrates is unconditional rather than conditional on having
survived to entry. That is usually what is wanted -- a lifetime or restricted
QALY is naturally measured from the time origin -- but it is a different
quantity from the one the likelihood conditioned on, so it is worth being
deliberate about. {cmd:jmqaly} does not currently offer the entry-conditional
version.

{marker examples}{...}
{title:Examples}

{pstd}Fit a joint model: a flexible parametric survival submodel sharing a
random intercept with a linear mixed model for quality of life{p_end}
{phang2}{cmd:. merlin (stime trt M1[id], family(rp, df(1) failure(died))) ///}{p_end}
{phang2}{cmd:          (qol time time#trt M1[id]@1, family(gaussian) timevar(time))}{p_end}

{pstd}A grid of horizons to predict over{p_end}
{phang2}{cmd:. range tvar 0 10 100}{p_end}

{pstd}Restricted QALYs under treatment{p_end}
{phang2}{cmd:. jmqaly qaly, at1(trt 1) timevar(tvar)}{p_end}

{pstd}Both arms, the difference, and confidence intervals for each{p_end}
{phang2}{cmd:. jmqaly dq, at1(trt 1) at2(trt 0) timevar(tvar) ci}{p_end}

{pstd}Discounted at 3.5% a year{p_end}
{phang2}{cmd:. jmqaly ddq, at1(trt 1) at2(trt 0) timevar(tvar) drate(0.035)}{p_end}

{pstd}Quality of life modelled on the logit scale, rescaled to the EQ-5D range{p_end}
{phang2}{cmd:. jmqaly eq, at1(trt 1) timevar(tvar) transform tomin(-0.594) tomax(1)}{p_end}

{pstd}The survival and RMST components of the same integrand{p_end}
{phang2}{cmd:. jmqaly s, at1(trt 1) timevar(tvar) survival}{p_end}
{phang2}{cmd:. jmqaly r, at1(trt 1) timevar(tvar) rmst}{p_end}

{pstd}With three submodels, name the one holding quality of life{p_end}
{phang2}{cmd:. jmqaly q, at1(trt 1) timevar(tvar) qolmodel(2)}{p_end}

{pstd}The lifetime QALY under treatment, bounding the utility scale so the
trajectory cannot be extrapolated below zero{p_end}
{phang2}{cmd:. jmqaly lt, lifetime at1(trt 1) drate(0.035) ///}{p_end}
{phang2}{cmd:          transform frommin(-3) frommax(3) tomin(0) tomax(1)}{p_end}
{phang2}{cmd:. display r(horizon)}{p_end}
{phang2}{cmd:. display r(remaining)}{p_end}

{pstd}The horizon it resolved and the contribution it estimates is still
missing are echoed, and returned in {cmd:r()}, because both belong with the
number in anything that reports it{p_end}

{marker references}{...}
{title:References}

{phang}
Crowther, M. J., A. Gasparini, S. Ekberg, F. Felizzi, E. Gallagher, and
N. Paracha. 2026. A Framework for the Estimation of Quality-Adjusted
Life-Years Using Joint Models of Longitudinal and Survival Data.
{it:Medical Decision Making}.
{browse "https://doi.org/10.1177/0272989X261477224"}

{phang}
Crowther, M. J. 2020. merlin - a unified framework for data analysis and
methods development in Stata. {it:Stata Journal} 20(4): 763-784.
{browse "https://journals.sagepub.com/doi/pdf/10.1177/1536867X20976311"}

{marker author}{...}
{title:Author}

{pstd}Michael J. Crowther and Alessandro Gasparini{break}
Red Door Analytics AB, Stockholm, Sweden{break}
{browse "mailto:michael.crowther@reddooranalytics.se":michael.crowther@reddooranalytics.se}
