// bivariate joint model: quality of life + survival
//
// run from the root of the jmqaly repo:
// 	do testing/testing_jmqaly.do

adopath ++ "`c(pwd)'/ado"
clear all
qui do "`c(pwd)'/mata/jmqaly.mata"

set seed 7254
set obs 300
gen id  = _n
gen u0  = rnormal(0,0.1)
gen u1  = rnormal(0,0.01)
gen trt = runiform()>0.5

survsim stime died, hazard(0.1:*1.2:*{t}:^0.2 :* 			///
        exp(-0.2 :* (0.7 :+ u0 :+ (-0.05:+u1 :+ 0.01:*trt):*{t})))	///
        covariates(trt -0.5) 						///
        maxt(10)							//

expand 10
bys id : gen time = _n-1
drop if time>stime

gen xb = 0.7 + u0 + (-0.05+u1 + 0.01*trt) * time
gen y  = rnormal(xb,0.5)

bys id (time) : replace stime 	= . if _n>1
bys id (time) : replace died 	= . if _n>1

qui merlin 	(stime trt M1[id]					///
		, family(rp, df(1) failure(died)))	///
	(y 	rcs(time, df(3) orthog) 				///
		time#trt						///
		M1[id]@1						///
		, family(gaussian) timevar(time))			//

range tvar 0 30 100

// restricted QALYs, and the difference between treatment arms
jmqaly qaly, at1(trt 1) timevar(tvar)
jmqaly dqaly, at1(trt 1) at2(trt 0) timevar(tvar) ci drate(0.02)

// conditional on the random effects being zero
jmqaly qfo, at1(trt 1) timevar(tvar) fixedonly

// lifetime QALYs: the horizon is resolved rather than given, and reported with
// the share of the total still estimated to be missing beyond it.
//
// this one WARNS, and that is the point of running it here: the toy quality of
// life trajectory above starts at 0.7 and falls by 0.05 a year, so it crosses
// zero inside the horizon the search reaches. a crossing is not a refusal --
// with discounting the integral still converges -- so jmqaly names the crossing
// point and returns the value.
jmqaly lt, at1(trt 1) at2(trt 0) lifetime drate(0.035)
di as txt "horizon `=r(horizon)', remaining `=r(remaining)'"

// the survival and RMST components, which should agree with -merlin- itself
jmqaly surv, at1(trt 1) timevar(tvar) survival
jmqaly rmst, at1(trt 1) timevar(tvar) rmst
predict msurv, survival marginal at(trt 1) outcome(1) timevar(tvar) zeros
predict mrmst, rmst marginal at(trt 1) outcome(1) timevar(tvar) zeros

gen double d_surv = reldif(surv_at1,msurv)
gen double d_rmst = reldif(rmst_at1,mrmst)

su qaly* dqaly* qfo*
su d_surv d_rmst
