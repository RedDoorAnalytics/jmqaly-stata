# Changelog

## v1.6.3   17sep2026

### Changed

- Engine internals removed from the changelog as well as the help: the
  description of how the compiled library binds to `merlin`, and the references
  to `merlin` internal functions. The requirement and its diagnostic are
  unchanged.

## v1.6.2   17sep2026

### Fixed

- **Two formulae in `help jmqaly` did not render.** `{integral}`, `{sub:}` and
  `{sup:}` are not SMCL directives, so they passed through literally and the
  QALY definition read

      QALY(t) = {integral}{sub:0}{sup:t} Q(u) S(u) exp(-r u) du

  and the discrete discount factor read `(1+r){sup:-u}`. Both now use
  characters the help viewer renders: `∫₀ᵗ` and `(1+r)⁻ᵘ`.

- Engine internals removed from the help: the section describing how the
  compiled library binds to `merlin`'s internal structures, and the one
  remaining reference to a `merlin` internal function. The version requirement
  and its diagnostic are unchanged.

## v1.6.1   17sep2026

### Fixed

- **The paper is now cited under its published title**, in the README and the
  help file: *A Framework for the Estimation of Quality-Adjusted Life-Years
  Using Joint Models of Longitudinal and Survival Data*, confirmed against
  Crossref for doi:10.1177/0272989X261477224. A working title had been carried
  in its place.

## v1.6.0   17sep2026

### Refused

- **A quality of life submodel fitted without `timevar()` is now refused.**
  `jmqaly` evaluates the trajectory at Gauss-Legendre nodes of its own
  choosing, and `merlin` can only re-evaluate a submodel at a new time when it
  knows which variable carries time there — which is what `timevar()` records.
  Fitted without it, the trajectory does not move with the node and the QALY
  integrates the wrong thing, **without erroring**: `merlin` returns a number
  and `jmqaly` integrates it.

  Measured on the example this README carries, both runs completing normally:

  | | QALY at *t*=10 | treatment difference |
  |---|---|---|
  | with `timevar(time)` | 5.11576075 | 1.54196013 |
  | without | 5.29318812 | 1.27785953 |

  2.8% on the QALY and 17% on the difference. A refusal is the only safe
  behaviour.

- **The survival submodel is deliberately not subject to this.** `merlin` sets
  its time variable from the response whether or not one is passed, and passing
  one where there is no time-dependent effect only forces `merlin` to integrate
  the cumulative hazard numerically instead of using the closed form. The
  examples in the README and the help no longer pass `timevar()` there.

## v1.5.2   17sep2026

### Fixed

- **A `lifetime` run that ran out of horizon too early was told its integral
  diverges.** Forming the increment ratio the search relies on takes three
  evaluations — one to set the previous value, one to set the previous
  increment, one to divide. A `maxhorizon()` admitting fewer produces no ratio
  at all, and that absent ratio was being reported as a *bad* one: the refusal
  claimed the integral was not converging and advised `drate()`, which the
  caller may already have set.

  Measured on a model that settles cleanly at horizon 128 with `drate(0.035)`:
  with the largest observed event time at 8, `maxhorizon(20)` admits two
  evaluations and forms no ratio. Nothing about that model diverges. It now
  reports how many evaluations fitted, says three are needed before increments
  can be compared, and states that this says nothing about convergence.

  The divergence diagnosis is unchanged and still fires on the case it is for.

## v1.5.1   16sep2026

### Changed

- `help jmqaly` no longer directs the reader to set the random-effect
  quadrature with `intpoints()` on the `merlin` command. The count is still
  `merlin`'s rather than `jmqaly`'s, and the help still says so.

## v1.5.0   16sep2026

### Changed

- **Requires Stata 19.5.** `jmqaly` declared `version 18` until now. This is
  `jmqaly`'s own floor and not one inherited from its engine: the released
  `merlin` 2.5.0 declares `version 15.1`, so anyone reading the previous
  README — which said 19.5 "is what `merlin` requires" — was reading something
  that is not true of the `merlin` this package ships against.
- **Alessandro Gasparini is credited as an author**, alongside Michael J.
  Crowther, in the package file, the help file and the README.
- The README no longer carries the long explanation of why the `merlin`
  version is pinned. The requirement stays; `help jmqaly` keeps the reasoning
  under *Which merlin*, which is where a user hits the refusal.

## v1.4.0   15sep2026

### Requires merlin 2.5.0

- **`jmqaly` now names the `merlin` it needs, and refuses any other.** The
  compiled library works only with the `merlin` it was built against, so a
  mismatch is checked for and refused with a diagnostic naming both versions,
  rather than left to fail later.

  **2.5.0 is the released `merlin`** — what `net install` from the `merlin`
  repository gives you.

  The `jmqaly` *source* is not version-specific; only the distributed library
  is bound. Rebuilding against another `merlin` is supported.

### New

- **Lifetime QALYs: `lifetime`, with `tolerance()` and `maxhorizon()`.**
  An unrestricted horizon was supported once and removed in April 2025.
  `lifetime` does not integrate to infinity: it doubles the horizon until the
  accumulated QALY stops moving, integrates to that point, and **reports the
  horizon and the contribution it estimates is still missing**, in `r(horizon)`
  and `r(remaining)`. "Integrated to 232, with 3e-05 of the total remaining" is
  a statement a reader can check; "integrated to infinity" is not.

  The search works on the estimand rather than on survival. Writing *d_k* for
  the increase over the *k*th doubling, the ratio *rho* = *d_k*/*d_k-1*
  measures the decay and the remainder is estimated as *d_k rho*/(1-*rho*). For
  a tail behaving like *t^-a* this ratio tends to 2^(1-*a*), below 1 exactly
  when *a* > 1 — so the accept/refuse test is the convergence condition itself
  rather than a proxy for it. A survival threshold would have been a proxy, and
  a poor one: under a shared random effect the marginal survival is far
  heavier-tailed than any individual's, so it can sit at 1e-4 over a long
  stretch while contributing almost nothing.

  The integral does not always exist. A model with a surviving fraction leaves
  an integrand that never vanishes and is refused undiscounted, naming
  `drate()`. A `tolerance()` below what the quadrature can demonstrate is
  refused separately, naming the noise floor — raising `qintpoints()` lowers it
  but does not remove it. A utility trajectory crossing zero is a **warning,
  not a refusal**: with discounting the integral still converges, so the value
  is returned with the crossing point named.

  `lifetime` may not be combined with `timevar()` — neither is a default,
  because a QALY curve and a lifetime QALY are different estimands — nor with
  `survival`, which does not accumulate. `lifetime rmst` works.

  Covered by `cert/cert_qaly_lifetime.do`, whose spine is an exact identity:
  the lifetime QALY must equal the restricted QALY at the horizon it reports,
  for `at1`, `at2` and the difference alike.

### Fixed

- A lifetime refusal no longer leaves `jmqaly_*` globals set for the next call.
  These refusals fire before the command's own error wrapper, so they had been
  bypassing the cleanup.

## v1.3.0

### New

- **Any survival family merlin can predict from.** The survival function is
  resolved from the fitted family with `merlin_p_getpf()` rather than naming
  the Royston-Parmar functions, so `weibull`, `exponential`, `gompertz`,
  `mspline`, `loghazard`, `logchazard` and the rest work alongside `rp`.
- **`dmethod(continuous|discrete)`** selects the discount factor: `exp(-r u)`
  as before, or `(1+r)^-u` as several health-technology guidelines specify.
- **Non-gaussian quality of life submodels.** The trajectory is taken on its
  natural scale, with the submodel's inverse link applied, so `bernoulli`,
  `gamma`, `poisson` and `beta` work as well as `gaussian`.

### Changed

- **The default `qintpoints` is 30, was 50.** The outer Gauss-Legendre
  integral is converged well before 50: at 30 nodes the QALY differs from a
  200-node evaluation by about 2e-08 relative. A `ci` prediction is about 11%
  faster as a result, since every node is re-evaluated for each perturbed
  parameter. Results move in the eighth decimal; pass `qintpoints(50)` to
  reproduce earlier output exactly.
- `fixedonly` works. It previously exited with a conformability error.
- `transform` respects `drate()`, and validates its from/to ranges instead of
  returning missing when they are left at their defaults.
- `drate()` is a real, is checked for sign, and is refused with `survival` or
  `rmst` rather than ignored.
- The QALY and the RMST are zero at t=0 rather than missing, and `survival`
  is 1 there with degenerate confidence limits, matching merlin's own
  `predict`. What that rests on for the log-time bases is declared in
  `cert/TOLERATED.txt`.
- Globals passed to the mata are `jmqaly_` prefixed and are dropped even when
  the prediction errors or is interrupted.

### Refused

- Cross-submodel association structures — `EV[]`, `dEV[]`, `iEV[]`,
  `d2EV[]` and `XB[]`. jmqaly evaluates the survival submodel at quadrature
  nodes of its own choosing, and these make the hazard at each node depend on
  another submodel evaluated at the same point, which that evaluation does not
  carry.
- Ordinal quality of life submodels, by name, since merlin supplies no
  expected value for them.

### Packaging

- `build/buildmlib.do` compiles `mata/` into `ljmqaly.mlib`. The command
  previously required the mata to be run by hand.
- `ado/jmqaly.sthlp`, `jmqaly.pkg` and `stata.toc`.
