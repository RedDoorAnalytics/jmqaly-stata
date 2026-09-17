*! version 1.6.3 17sep2026 MJC & AG

/*
- postestimation command following a joint model estimated with merlin
- calculates restricted quality adjusted life years

	QALY(t) = int_0^t Q(u) S(u) D(u) du

  with D(u) the discount factor: exp(-drate*u) by default, or the discrete
  (1+drate)^(-u) that some health technology guidelines specify

  where Q(u) is the quality of life trajectory and S(u) the survival
  probability, both taken from the joint model in memory. The integral over
  time is Gauss-Legendre; the integral over the random effects is merlin's.

- the survival function is resolved from the fitted family, so any survival
  family merlin can predict from is supported
- the fit must be two-level with exactly one survival submodel
- EV[] and XB[] association structures are not supported
*/


//----------------------------------------------------------------------------//
// Drop every global the mata functions read. Called before ANY exit, because a
// lifetime refusal happens before the main capture-noisily wrapper and would
// otherwise leave jmqaly_* set for whatever runs next in the session.
//----------------------------------------------------------------------------//
program define jmqaly_cleanup
	capture macro drop jmqaly_stmodel jmqaly_xtmodel jmqaly_qintpoints
	capture macro drop jmqaly_drate jmqaly_ddisc
	capture macro drop jmqaly_transform
	capture macro drop jmqaly_frommin jmqaly_frommax jmqaly_tomin jmqaly_tomax
	capture drop __jmq_q
end

//----------------------------------------------------------------------------//
// One evaluation of the QALY at a fixed horizon, without CIs. Used by the
// lifetime horizon search, which needs the integrand's own behaviour rather
// than a proxy for it.
//----------------------------------------------------------------------------//
program define jmqaly_lt_value, rclass
	syntax , QFunc(string) HZ(real) PREDType(string) [ ATspec(string) ZEROs ]
	tempvar tv
	qui gen double `tv' = `hz'
	local atopt
	if `"`atspec'"'!="" local atopt at(`atspec')
	capture drop __jmq_q
	qui predict double __jmq_q, userfunction(`qfunc') timevar(`tv') ///
		`predtype' `zeros' `atopt'
	qui summarize __jmq_q, meanonly
	return scalar v = r(mean)
	capture drop __jmq_q
end

//----------------------------------------------------------------------------//
// Resolve the lifetime horizon on the ESTIMAND, not on survival.
//
// A survival tolerance is a proxy, and it comes apart from the thing we care
// about exactly where it matters: under a shared frailty the MARGINAL survival
// has a far heavier tail than any individual's, so S can sit at 1e-4 over a
// long stretch while contributing almost nothing to the QALY. Bounding the
// remaining CONTRIBUTION instead dissolves that, and needs no assumption that
// utility is bounded.
//
// Double the horizon until the accumulated QALY stops moving. With three
// successive values the increments d_k give a ratio rho = d_k/d_{k-1}; while
// rho < 1 the integrand is decreasing geometrically and the remainder beyond
// the current horizon is estimated by d_k * rho/(1-rho). Stop when that
// remainder is below tolerance() of the accumulated value, and REPORT it.
//
// rho >= 1 means the increments are not shrinking: either a surviving
// fraction, whose undiscounted integral genuinely diverges, or a horizon still
// far too short. Those need opposite advice, so the caller is told which.
//----------------------------------------------------------------------------//
program define jmqaly_horizon, rclass
	syntax , QFunc(string) T0(real) TOLer(real) MAXH(real) ///
		PREDType(string) [ ATspec(string) ZEROs ]

	local hz     = `t0'
	local dprev  = .
	local vprev  = .
	//HOW MANY EVALUATIONS ACTUALLY FITTED. Forming an increment RATIO needs
	//three: the first sets vprev, the second sets dprev, the third divides.
	//A maxhorizon() admitting fewer yields no ratio at all, which is a
	//different state from a ratio that is bad, and the caller has to be able
	//to tell them apart.
	local neval  = 0
	local v      = .
	local rho    = .
	local sawneg = .
	local lastneg = 0
	local nbad   = 0
	local bestinc = .
	local firstinc = .

	// Do NOT accept or refuse on a single rho. Three successive values put
	// the increments in their asymptotic regime only if the integrand is
	// already decaying; a utility trajectory that RISES before it falls
	// makes the increment grow first, so an early rho >= 1 says the
	// geometric regime has not started, not that the integral diverges.
	// Keep doubling and let the evidence accumulate; only the state at
	// maxhorizon is diagnostic. Same discipline as probing twice, one
	// level up.
	while `hz' <= `maxh' {
		jmqaly_lt_value , qfunc(`qfunc') hz(`hz') ///
			predtype(`predtype') atspec(`atspec') `zeros'
		local v = r(v)
		local ++neval
		if missing(`v') {
			return scalar hz  = .
			return scalar bad = 1
			return scalar v   = .
			exit
		}
		if !missing(`vprev') {
			local d = `v' - `vprev'
			// The ratio is only meaningful while the increments are.
			// Once they reach the numerical noise floor the QALY has
			// stopped moving but rho becomes garbage -- traced at
			// 1.49 then 6.77 on a model whose value had changed by
			// 1e-4 relative. Chasing rho past that point refuses a
			// converged integral, so test the increment's own size
			// first and prefer the geometric estimate only while it
			// can be formed.
			local relinc = cond(`v'!=0, abs(`d'/`v'), .)
			// smallest relative increment achieved anywhere in the
			// search: this is what separates 'the model diverges'
			// from 'you asked for a tolerance below the quadrature
			// noise floor', which need opposite advice.
			if !missing(`relinc') & missing(`firstinc') {
				local firstinc = `relinc'
			}
			if !missing(`relinc') {
				if missing(`bestinc') {
					local bestinc = `relinc'
				}
				else if `relinc' < `bestinc' {
					local bestinc = `relinc'
				}
			}
			if !missing(`dprev') & `dprev'!=0 {
				local rho = `d'/`dprev'
				if `rho' >= 0 & `rho' < 1 {
					local rem = `d'*`rho'/(1-`rho')
					local relrem = cond(`v'!=0, abs(`rem'/`v'), .)
					if `relrem' < `toler' {
						local lastneg = 0
						local nbad = 0
						return scalar hz      = `hz'
						return scalar value   = `v'
						return scalar relrem  = `relrem'
						return scalar rho     = `rho'
						return scalar sawneg  = `sawneg'
						exit
					}
				}
				else if `relinc' < `toler' {
					// noise floor: rho unusable, but the integral
					// has demonstrably stopped moving
					return scalar hz      = `hz'
					return scalar value   = `v'
					return scalar relrem  = `relinc'
					return scalar rho     = .
					return scalar sawneg  = `sawneg'
					exit
				}
				if `d' < 0 {
					// integrand has changed sign in this window
					if missing(`sawneg') local sawneg = `hz'
					local lastneg = 1
					local nbad = `nbad' + 1
				}
				else if `rho' >= 1 {
					local lastneg = 0
					local nbad = `nbad' + 1
				}
				else {
					local lastneg = 0
					local nbad = 0
				}
			}
			local dprev = `d'
		}
		local vprev = `v'
		local hz = 2*`hz'
	}
	// exhausted the horizon: the state HERE is what diagnoses it
	return scalar hz      = .
	return scalar neval   = `neval'
	return scalar rho     = `rho'
	return scalar v       = `v'
	return scalar nbad    = `nbad'
	return scalar lastneg = `lastneg'
	return scalar sawneg  = `sawneg'
	return scalar bestinc  = `bestinc'
	return scalar firstinc = `firstinc'
end


program define jmqaly, rclass
	version 19.5
	syntax [anything(id="newvarname")] 				///
			, 	[					///
				TIMEVar(varname)			///
				LIFEtime				///
				TOLerance(real 1e-4)			///
				MAXHorizon(real 0)			///
				AT1(string)				///
				AT2(string)				///
				QOLModel(numlist int min=1 >0)	 	///
				QINTPOINTS(numlist int min=1 >0) 	///
				FIXEDONLY				///
				STANDardise				///
				SURVival				///
				RMST 					///
				CI					///
				ITERate(real 100)			///
									///
				frommin(real 0)				///
				frommax(real 1)				///
				tomin(real 0)				///
				tomax(real 1)				///
				TRANSFORM				///
									///
				DRATE(real 0)				///
				DMEthod(string)				///
			]

	// exactly one of timevar() and lifetime. Neither is a default: a QALY
	// curve and a lifetime QALY are different estimands and the caller
	// should say which is wanted.
	if "`timevar'"=="" & "`lifetime'"=="" {
		di as error "specify timevar() for a QALY over a time grid, " ///
			"or lifetime for the whole horizon"
		exit 198
	}
	if "`lifetime'"!="" & "`survival'"!="" {
		di as error "lifetime may not be combined with survival: a " ///
			"survival probability does not accumulate,"
		di as error "so there is no remaining contribution to bound and " ///
			"no horizon to resolve. Use"
		di as error "timevar() for survival at chosen times, or lifetime " ///
			"with rmst for time alive."
		exit 198
	}
	if "`timevar'"!="" & "`lifetime'"!="" {
		di as error "timevar() and lifetime are mutually exclusive"
		exit 198
	}
	if "`lifetime'"=="" & (`maxhorizon'!=0 | `tolerance'!=1e-4) {
		di as error "tolerance() and maxhorizon() require lifetime"
		exit 198
	}
	if `tolerance'<=0 | `tolerance'>=1 {
		di as error "tolerance() must be in (0,1)"
		exit 198
	}

	// error checks
	if "`e(cmd)'"!="merlin" {
		di as error "Joint model not in memory"
		exit 198
	}

	// WHICH merlin THIS COPY OF jmqaly WAS COMPILED AGAINST.
	//
	// The compiled ljmqaly.mlib is bound to the merlin it was built against,
	// not to a stable interface, and the binding changes between merlin
	// releases. Run the library against a different merlin and it fails inside
	// merlin, with an error naming neither package, neither version, and no
	// cause. Refusing here costs one string comparison and turns that into a
	// sentence a user can act on. The same jmqaly SOURCE works against every
	// merlin tried; it is the compiled library that is bound, so rebuilding is
	// a real remedy and the message says so.
	//
	// The supported version is declared here and asserted by
	// cert_qaly_engine.do against the merlin a cert run actually resolved, so
	// this line and the artefact cannot drift apart silently.
	local JMQALY_MERLIN 2.5.0
	local _mver "unresolved"
	capture noisily {
		findfile merlin.ado
		tempname _mfh
		file open `_mfh' using `"`r(fn)'"', read
		file read `_mfh' _mline
		file close `_mfh'
		local _p = strpos(`"`_mline'"', "version")
		if `_p' {
			local _rest = substr(`"`_mline'"', `_p'+8, .)
			local _mver : word 1 of `_rest'
		}
	}
	if "`_mver'" != "`JMQALY_MERLIN'" {
		di as error "{p}this jmqaly was built against merlin `JMQALY_MERLIN', " ///
			"but the merlin in use reports `_mver'.{p_end}"
		di as error "{p}jmqaly's compiled library is bound to merlin's " ///
			"internal structures, which change between merlin releases, " ///
			"so it runs only with the merlin it was built against. With " ///
			"any other it fails inside merlin itself, with an error " ///
			"naming neither package.{p_end}"
		di as error "{p}Either install merlin `JMQALY_MERLIN', or rebuild " ///
			"jmqaly from source against the merlin you have -- " ///
			"build/buildmlib.do in the jmqaly repository.{p_end}"
		exit 198
	}
	if `e(Nmodels)'<2 {
		di as error "At least two submodels must be present"
		exit 198
	}
	if `e(Nlevels)'!=2 {
		di as error "Only two-level joint models are supported"
		exit 198
	}
	local Nstmodels = 0
	forvalues i=1/`e(Nmodels)' {
		if "`e(failure`i')'"!="" {
			local Nstmodels = `Nstmodels' + 1
			local jmqaly_stmodel = `i'
		}
	}
	if !`Nstmodels' {
		di as error "Survival submodel not found"
		exit 198
	}
	if `Nstmodels'>1 {
		di as error "Only one survival submodel can be specified"
		exit 198
	}

	// The survival submodel is evaluated at quadrature nodes of jmqaly's own
	// choosing. A cross-submodel association makes the hazard at u depend on
	// another submodel evaluated at u, which is not carried through that
	// evaluation, so those fits are refused rather than answered wrongly.
	// "EV[" matches EV[], dEV[], iEV[] and d2EV[] alike; XB[] is the same
	// mechanism reading the linear predictor instead of the expected value.
	forvalues i=1/`e(Nmodels)' {
		local cl "`e(cmplabels`i')'"
		if strpos("`cl'","EV[") | strpos("`cl'","XB[") {
			di as error ///
	"EV[] and XB[] association structures are not supported"
			di as error ///
	"  submodel `i' uses one; jmqaly needs an association that does not"
			di as error ///
	"  depend on another submodel evaluated at the same time point"
			exit 198
		}
	}

	if `e(Nmodels)'>2 & "`qolmodel'"=="" {
		di as error "qolmodel() required"
		exit 198
	}

	if `e(Nmodels)'==2 & "`qolmodel'"=="" {
		local jmqaly_xtmodel = 1
		if `jmqaly_stmodel'==1 {
			local jmqaly_xtmodel = 2
		}
	}
	else {
		local jmqaly_xtmodel = `qolmodel'
	}

	// THE QUALITY OF LIFE SUBMODEL MUST CARRY timevar().
	//
	// jmqaly evaluates the trajectory at Gauss-Legendre nodes of its OWN
	// choosing, not at the times in the data. merlin can only re-evaluate a
	// submodel at a new time when it knows which variable carries time there,
	// and that is what timevar() records -- e(timevar#) is empty without it.
	// Fitted without it, the trajectory does not move with the node and the
	// QALY integrates the wrong thing. It does not error: merlin returns a
	// number, jmqaly integrates it, and the answer is quietly different.
	//
	// MEASURED on the example this package's own README carries, both runs
	// completing with rc 0:
	//     with    timevar(time)   QALY(10) 5.11576075   difference 1.54196013
	//     without timevar(time)   QALY(10) 5.29318812   difference 1.27785953
	// 2.8% on the QALY and 17% on the treatment difference.
	//
	// The SURVIVAL submodel is a different case and deliberately not checked:
	// merlin sets its timevar from the response, so e(timevar#) is populated
	// whether or not the user passed one -- and passing one there, with no
	// time-dependent effect, only forces merlin to integrate the cumulative
	// hazard numerically instead of using the closed form.
	//
	// Scoped to the QALY. survival and rmst never touch the quality of life
	// submodel, so refusing them for its timevar would be refusing a fit that
	// is perfectly adequate for what was asked.
	if "`survival'"=="" & "`rmst'"=="" {
		if "`e(timevar`jmqaly_xtmodel')'"=="" {
			di as error ///
	"the quality of life submodel was fitted without timevar()"
			di as error ///
	"  submodel `jmqaly_xtmodel' has no time variable recorded, so merlin cannot"
			di as error ///
	"  re-evaluate the trajectory at the quadrature nodes the QALY integrates"
			di as error ///
	"  over, and the result would be silently wrong rather than an error."
			di as error ///
	"  Refit with timevar() on that submodel, naming the variable that carries"
			di as error ///
	"  time in it. The survival submodel does not need one."
			exit 198
		}
	}

	// merlin supplies an expected value only for the families listed in
	// merlin_setup_EV_p(); for any other the pointer it hands back is not a
	// function, and merlin_util_expval() fails deep in mata with "matrix
	// found where function required". Refuse here instead, by name. This
	// list mirrors merlin's and has to grow with it.
	local evfams gaussian bernoulli gamma poisson beta null
	local qolfam "`e(family`jmqaly_xtmodel')'"
	if !`: list qolfam in evfams' {
		di as error ///
	"family(`qolfam') quality of life submodels are not supported"
		di as error ///
	"  merlin has no expected value for family(`qolfam'); jmqaly needs one"
		di as error ///
	"  supported: `evfams'"
		exit 198
	}

	if "`survival'"!="" & "`rmst'"!="" {
		di as error "only one of 'survival' and 'rmst' can be selected"
		exit 198
	}

	if "`anything'"=="" {
		di as error "newvarname required"
		exit 100
	}

	// drate()/transform() only apply to the QALY calculation
	if `drate'<0 {
		di as error "drate() must be non-negative"
		exit 198
	}
	if "`dmethod'"=="" {
		local dmethod continuous
	}
	if !inlist("`dmethod'","continuous","discrete") {
		di as error ///
	"dmethod() must be continuous or discrete"
		exit 198
	}
	if `drate' & ("`survival'"!="" | "`rmst'"!="") {
		di as error "drate() may not be combined with 'survival' or 'rmst'"
		exit 198
	}
	if "`transform'"!="" & ("`survival'"!="" | "`rmst'"!="") {
		di as error ///
		"'transform' may not be combined with 'survival' or 'rmst'"
		exit 198
	}
	if "`transform'"!="" {
		if `frommax'<=`frommin' {
			di as error "frommax() must be greater than frommin()"
			exit 198
		}
		if `tomax'<=`tomin' {
			di as error "tomax() must be greater than tomin()"
			exit 198
		}
	}

	confirm new var `anything'_at1
	if "`ci'"!="" {
		confirm new var `anything'_at1_lci
		confirm new var `anything'_at1_uci
	}
	if "`at2'"!="" {
		confirm new var `anything'_at2
		if "`ci'"!="" {
			confirm new var `anything'_at2_lci
			confirm new var `anything'_at2_uci
		}
		confirm new var `anything'_diff
		if "`ci'"!="" {
			confirm new var `anything'_diff_lci
			confirm new var `anything'_diff_uci
		}
	}

	// pass settings through to the mata functions
	global jmqaly_stmodel = `jmqaly_stmodel'
	global jmqaly_xtmodel = `jmqaly_xtmodel'

	if "`qintpoints'"=="" {
		local qintpoints = 30
	}
	global jmqaly_qintpoints = `qintpoints'

	if `drate' {
		global jmqaly_drate = `drate'
		if "`dmethod'"=="discrete" {
			global jmqaly_ddisc = 1
		}
	}
	if "`transform'"!="" {
		global jmqaly_transform = "transform"
		global jmqaly_frommin = `frommin'
		global jmqaly_frommax = `frommax'
		global jmqaly_tomin   = `tomin'
		global jmqaly_tomax   = `tomax'
	}

	// get stuff
	local zeros zeros
	if "`standardise'"!="" {
		local zeros
	}

	local predtype "marginal"
	if "`fixedonly'"!="" {
		local predtype "fixedonly"
	}

	// label according to prediction type
	local lbl = "QALY"
	local qfunc jmqaly_qaly
	if "`survival'"!="" {
		local lbl = "S"
		local qfunc jmqaly_surv
	}
	if "`rmst'"!="" {
		local lbl = "RMST"
		local qfunc jmqaly_rmst
	}

	// ---- lifetime: resolve a horizon and integrate to it -------------//
	// Integrating to the point where the cohort is effectively dead, and
	// REPORTING that point, rather than claiming an infinite integral. The
	// horizon is common to both arms -- using arm-specific horizons would
	// mix a horizon effect into the contrast.
	if "`lifetime'"!="" {
		local rv : word 1 of `e(response`jmqaly_stmodel')'
		qui summarize `rv', meanonly
		local tmax = r(max)
		if `maxhorizon'<=0 {
			local maxhorizon = 500*`tmax'
		}
		local hz = 0
		local relrem = 0
		local negfrom = .
		foreach a in 1 2 {
			if `a'==2 & "`at2'"=="" continue
			local thisat "`at1'"
			if `a'==2 local thisat "`at2'"
			jmqaly_horizon , qfunc(`qfunc') t0(`tmax')            ///
				toler(`tolerance') maxh(`maxhorizon')          ///
				predtype(`predtype') atspec(`thisat') `zeros'
			if missing(r(hz)) {
				// An evaluation that came back MISSING is a failure to
				// evaluate, not a statement about convergence. Saying
				// "your integral diverges" here would be a confident
				// diagnosis of something never measured.
				if r(bad)==1 {
					di as error "the QALY could not be evaluated at a " ///
						"horizon the search needed."
					di as error "This is an evaluation failure rather " ///
						"than a statement about convergence:"
					di as error "the model returned missing where it was " ///
						"asked to predict. Restrict the"
					di as error "horizon with timevar(), or lower " ///
						"maxhorizon()."
					jmqaly_cleanup
					exit 459
				}
				local rho = r(rho)
				local bi  = r(bestinc)
				local fi  = r(firstinc)
				//read BEFORE anything else touches r(): neval is the count
				//the too-few-evaluations message reports, and the message is
				//useless without it
				local nev = r(neval)
				// ORDER MATTERS. The noise floor produces NEGATIVE
				// increments on an integrand that never changes sign,
				// so testing the sign first advises transform() for a
				// quadrature artefact. A negative increment only means
				// a sign-changing integrand while the increments are
				// still MEANINGFUL, which is the same precondition rho
				// needs. So: settled-or-noise first, sign second.
				//
				// Assumes increments that fall and then stop falling.
				// An integrand whose increments fall and then grow
				// without bound -- unbounded utility outrunning the
				// survival decay -- would read as "settled". transform()
				// bounds the utility and rules that out; documented
				// rather than defended against.
				if !missing(`bi') & !missing(`fi') & `bi' < 0.01*`fi' {
					di as error "tolerance(" %9.0g `tolerance' ") is below " ///
						"what this quadrature can demonstrate."
					di as error "The integral settled to within " %9.2e `bi' ///
						" of its value and then stopped"
					di as error "improving - further doublings move it only " ///
						"by numerical noise."
					di as error "Loosen tolerance(), or raise qintpoints() " ///
						"to LOWER the noise floor (it does not remove it)."
					jmqaly_cleanup
				exit 459
				}
				di as error "the lifetime integral has not settled by " ///
					"maxhorizon(" %9.0g `maxhorizon' ")."
				//THREE STATES, NOT TWO. An ABSENT increment ratio was
				//being reported as a BAD one: this read
				//`if r(nbad) > 0 | missing(rho)', so a maxhorizon() that
				//admitted fewer than the three evaluations a ratio needs
				//fell into the divergence branch and told the user their
				//integral was not converging -- advising drate(), which
				//they may already have set. Measured on a model that
				//settles at horizon 128 with drate(0.035): maxhorizon(20)
				//from tmax 8 admits two evaluations, forms no ratio, and
				//nothing whatever about it diverges.
				//
				//Order matters: absent is tested BEFORE bad, because a
				//missing rho makes every comparison against it false.
				if missing(`rho') {
					di as error "maxhorizon(" %9.0g `maxhorizon' ") admitted " ///
						"only " %3.0f `nev' " evaluation(s) from t = " ///
						%9.0g `tmax' ","
					di as error "and three are needed before the increments " ///
						"can be compared at all."
					di as error "This says nothing about whether the integral " ///
						"converges. Raise maxhorizon()."
				}
				else if r(nbad) > 0 {
					di as error "Successive doublings of the horizon are " ///
						"adding as much as the previous one"
					di as error "(ratio " %6.3f `rho' "), so the integral is " ///
						"NOT converging. A model with a"
					di as error "surviving fraction has a lifetime QALY that " ///
						"diverges without discounting:"
					di as error "use drate(), or restrict the horizon with " ///
						"timevar()."
				}
				else {
					di as error "It is converging (increment ratio " %6.3f `rho' ///
						") but has not yet reached"
					di as error "tolerance(" %9.0g `tolerance' "). Raise " ///
						"maxhorizon(), or loosen tolerance()."
				}
				jmqaly_cleanup
				exit 459
			}
			if r(hz) > `hz' {
				local hz = r(hz)
			}
			if !missing(r(sawneg)) & missing(`negfrom') {
				local negfrom = r(sawneg)
			}
			if r(relrem) > `relrem' {
				local relrem = r(relrem)
			}
		}
		// A utility trajectory that goes negative does NOT make the
		// integral diverge -- with discounting v*Q*S still vanishes,
		// and the lifetime QALY converges to a value that includes the
		// negative stretch. Refusing it was too strong; accepting it
		// silently is too weak. Warn, and name the horizon, which is
		// what a reader needs to judge whether the extrapolation is
		// meaningful. Genuine divergence is caught separately by the
		// increments failing to shrink.
		if !missing(`negfrom') {
			di as txt "{p}note: the utility trajectory is negative " ///
				"beyond t = " %9.0g `negfrom' ", so the lifetime "  ///
				"QALY includes negative utility from there. "       ///
				"Bound the utility scale with transform if that "   ///
				"is not intended.{p_end}"
		}
		tempvar _lthz
		qui gen double `_lthz' = `hz'
		local timevar `_lthz'
		di as txt "lifetime horizon: " as result %9.0g `hz'          ///
			as txt "   estimated remaining contribution: "        ///
			as result %6.3g 100*`relrem' as txt "% of the total"
	}

	local predopts	userfunction(`qfunc') 	///
			`zeros' 		///
			`standardise'		///
			timevar(`timevar')	///
			`predtype'
	if "`at1'"!="" {
		local at1opt at(`at1')
	}
	if "`at2'"!="" {
		local at2opt at(`at2')
	}

	// prediction
	capture noisily {

		quietly {
			if "`ci'"!="" {
				local cis ci(`anything'_at1_lci `anything'_at1_uci)
			}
			predictnl double `anything'_at1 = 	///
				predict(`predopts' `at1opt')	///
				, `cis' iterate(`iterate')
			label var `anything'_at1 "`lbl'(t | at1)"
			if "`ci'"!="" {
				label var `anything'_at1_lci ///
					"Lower 95% CI of `lbl'(t | at1)"
				label var `anything'_at1_uci ///
					"Upper 95% CI of `lbl'(t | at1)"
			}
		}

		if "`at2'"!="" {
			quietly {
				if "`ci'"!="" {
					local cis 			///
					ci(`anything'_at2_lci `anything'_at2_uci)
				}
				predictnl double `anything'_at2 = 	///
					predict(`predopts' `at2opt')	///
					, `cis' iterate(`iterate')
				label var `anything'_at2 "`lbl'(t | at2)"
				if "`ci'"!="" {
					label var `anything'_at2_lci ///
					"Lower 95% CI of `lbl'(t | at2)"
					label var `anything'_at2_uci ///
					"Upper 95% CI of `lbl'(t | at2)"
				}
			}

			quietly {
				if "`ci'"!="" {
					local diffcis 			///
					ci(`anything'_diff_lci `anything'_diff_uci)
				}
				predictnl double `anything'_diff = 	///
					predict(`predopts' `at1opt')	///
					- 				///
					predict(`predopts' `at2opt')	///
					, `diffcis' iterate(`iterate')
				label var `anything'_diff 	///
					"`lbl'(t | at1) - `lbl'(t | at2)"
				if "`ci'"!="" {
					label var `anything'_diff_lci ///
				"Lower 95% CI of `lbl'(t | at1) - `lbl'(t | at2)"
					label var `anything'_diff_uci ///
				"Upper 95% CI of `lbl'(t | at1) - `lbl'(t | at2)"
				}
			}
		}

		// S(0)=1. merlin's predict substitutes this after the prediction
		// is formed, in its ado layer, which a user function reached
		// through predictnl never sees -- so without this, jmqaly and
		// predict disagree at the origin for the families whose
		// cumulative hazard is a spline in log(t). Mirrors merlin as of
		// 5c8ef86, including the degenerate confidence limits. What it
		// rests on -- the sign of the fitted spline's lower-tail slope,
		// which is not checkable at runtime outside a baseline-only rp
		// model -- is stated in the help file under Remarks.
		if "`survival'"!="" {
			quietly {
				replace `anything'_at1 = 1 if `timevar'==0
				if "`ci'"!="" {
					replace `anything'_at1_lci = 1 ///
						if `timevar'==0
					replace `anything'_at1_uci = 1 ///
						if `timevar'==0
				}
				if "`at2'"!="" {
					replace `anything'_at2 = 1 ///
						if `timevar'==0
					replace `anything'_diff = 0 ///
						if `timevar'==0
					if "`ci'"!="" {
						replace `anything'_at2_lci = 1 ///
							if `timevar'==0
						replace `anything'_at2_uci = 1 ///
							if `timevar'==0
						replace `anything'_diff_lci = 0 ///
							if `timevar'==0
						replace `anything'_diff_uci = 0 ///
							if `timevar'==0
					}
				}
			}
		}

	}
	local rc = _rc

	// the horizon and the survival at it go into r(): these numbers are
	// carried into a submission and should not have to be read off the
	// console. Same reasoning as echoing the discount rate.
	if "`lifetime'"!="" {
		return scalar horizon   = `hz'
		return scalar remaining = `relrem'
	}

	capture macro drop jmqaly_stmodel jmqaly_xtmodel jmqaly_qintpoints
	capture macro drop jmqaly_drate jmqaly_ddisc
	capture macro drop jmqaly_transform
	capture macro drop jmqaly_frommin jmqaly_frommax jmqaly_tomin jmqaly_tomax

	if `rc' {
		exit `rc'
	}

end
