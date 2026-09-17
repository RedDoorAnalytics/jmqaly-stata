# jmqaly

Quality-adjusted life years after a joint longitudinal-survival model fitted
with [`merlin`](https://github.com/RedDoorAnalytics/merlin) in Stata.

`jmqaly` is the reference implementation of the estimators in:

> Crowther, M. J., A. Gasparini, S. Ekberg, F. Felizzi, E. Gallagher, and
> N. Paracha (2026). A Framework for the Estimation of Quality-Adjusted
> Life-Years Using Joint Models of Longitudinal and Survival Data.
> *Medical Decision Making*.
> [doi:10.1177/0272989X261477224](https://doi.org/10.1177/0272989X261477224)

## What it computes

For a horizon *t*, the restricted QALY

    QALY(t) = ∫₀ᵗ Q(u) S(u) exp(-r u) du

where `Q(u)` is the mean quality of life trajectory, `S(u)` the survival
probability, and `r` an optional discount rate. The integral over time is
Gauss-Legendre; the integral over the random effects is merlin's.

Q and S are multiplied at each quadrature node of the random effect
distribution and the product is then integrated over that distribution, so the
estimand is `E[Q(u)S(u)]` — it carries `Cov(Q,S)`, the association between
quality of life and survival that the joint model estimates.

With `lifetime` in place of `timevar()`, the horizon is resolved rather than
given: it doubles until the accumulated QALY stops moving, and the result is
reported with the horizon reached and the share of the total still estimated to
be missing beyond it. "Integrated to 232, with an estimated 3e-05 of the total
remaining" is a statement a reader can check; "integrated to infinity" is not.

## Requirements

- **StataNow 19.5.**
- **[`merlin`](https://github.com/RedDoorAnalytics/merlin) 2.5.0**, and that
  version specifically. It is the released `merlin`, so it is what
  `net install` from the `merlin` repository gives you. `jmqaly` checks for it
  and refuses to run beside any other.
- `survsim`, for the testing scripts.

## Install

```stata
net install jmqaly, from("https://raw.githubusercontent.com/RedDoorAnalytics/jmqaly-stata/main/")
```

To run this development tree instead, put it on the adopath and index the
compiled library — prefer a clean `adopath ++` over a stale `net install` copy,
which can desync the ado from `ljmqaly.mlib`:

```stata
local JMQALY /path/to/jmqaly
adopath ++ "`JMQALY'/ado"
adopath ++ "`JMQALY'"
mata: mata mlib index
```

## Use

```stata
// a joint model: flexible parametric survival sharing a random intercept
// with a linear mixed model for quality of life
merlin (stime trt M1[id], family(rp, df(1) failure(died))) ///
       (qol time time#trt M1[id]@1, family(gaussian) timevar(time))

range tvar 0 10 100

// restricted QALYs by arm, the difference, and confidence intervals
jmqaly dq, at1(trt 1) at2(trt 0) timevar(tvar) ci

// discounted at 3.5% a year, continuously or on the discrete convention
jmqaly ddq, at1(trt 1) at2(trt 0) timevar(tvar) drate(0.035)
jmqaly dsq, at1(trt 1) at2(trt 0) timevar(tvar) drate(0.035) dmethod(discrete)

// quality of life modelled on the logit scale, rescaled to the EQ-5D range
jmqaly eq, at1(trt 1) timevar(tvar) transform tomin(-0.594) tomax(1)

// lifetime QALYs: the horizon is found rather than given, and reported
jmqaly lt, at1(trt 1) at2(trt 0) lifetime drate(0.035)
display "integrated to " r(horizon) ", with " r(remaining) " of the total remaining"
```

See `help jmqaly` for the full syntax.

## Scope

- The survival function is resolved from the fitted family, so any survival
  family merlin can predict from is supported.
- Cross-submodel associations are refused — `EV[]`, `dEV[]`, `iEV[]`,
  `d2EV[]` and `XB[]`; share the random effects between submodels instead.
- `survival` is 1 at t=0, matching merlin's own `predict`. For `rp` and
  `logchazard` that rests on the sign of the fitted spline's lower-tail
  slope, which holds under merlin's default knots. The QALY and RMST are zero
  at t=0 for every family, which is a theorem.
- Under `ltruncated()`, the QALY still accrues from time zero and the
  survival integrated is unconditional — not conditional on surviving to
  entry. A different quantity from the one the likelihood conditioned on.
- Quality of life may be `gaussian`, `bernoulli`, `gamma`, `poisson` or
  `beta`. merlin supplies no expected value for its ordinal families, so those
  are refused.
- The fit must be two-level, with exactly one survival submodel.
- A QALY is either restricted to a horizon you give in `timevar()`, or
  `lifetime`, which resolves one and reports it. Neither is a default: a QALY
  curve and a lifetime QALY are different estimands. A `lifetime` horizon
  extrapolates the fitted survival well beyond the data with no
  population-rate anchor, so the result is a model quantity rather than an
  all-cause one.
- `standardise` needs merlin support for standardisation after a model with
  random effects.

## Repository layout

    ado/          jmqaly.ado and its help file
    mata/         the mata source
    build/        buildmlib.do, which compiles mata/ into ljmqaly.mlib
    testing/      an end-to-end example, run from the repo root

## Building

```sh
/path/to/stata-mp -b do build/buildmlib.do
```

`ljmqaly.mlib` is committed, so a checkout is directly installable; rebuild it
after any change under `mata/`.

The rebuild binds the library to whichever `merlin` is on the ado-path when it
runs — `ljmqaly.mlib` is compiled against `merlin`'s internal structures, so it
works only with the `merlin` it was built against. If you rebuild against a
`merlin` other than 2.5.0, change the version `ado/jmqaly.ado` declares in
`JMQALY_MERLIN` to match, or the runtime check will refuse the `merlin` you
just built for.

## Certification

`jmqaly` has a certification suite of nine scripts, with pass/fail criteria
recorded and dated before execution, a register of the checks it deliberately
does not perform, and a requirement-to-cert traceability matrix. The flagship
cert validates the marginal QALY against a Monte Carlo truth computed from the
data generating process, and additionally requires a deliberately wrong
implementation — one that marginalises quality of life and survival separately
— to fail the same criterion. A suite that passes both certifies nothing.

The suite is evidence offered under a qualification pack rather than part of
the distribution, so it is not in this repository.

## Author and licence

Michael J. Crowther and Alessandro Gasparini, Red Door Analytics AB, Stockholm.

Licence: **GPL-3.0-or-later** — see [`LICENSE`](LICENSE).
Copyright (C) 2024-2026 Red Door Analytics.
