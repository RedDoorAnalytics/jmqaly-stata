version 19.5

local RS 	real scalar
local RM 	real matrix
local RC 	real colvector
local PS 	pointer scalar
local JM 	struct merlin_struct scalar

mata:

/*
Each function below returns predictions CONDITIONAL on the random effects,
i.e. an N x ndim matrix with one column per quadrature node. -merlin-
integrates over the random effect distribution itself (predict, marginal);
under fixedonly there is no node dimension and nothing to integrate.

The estimand is

	QALY(t) = int_0^t E_b[ Q(u,b) S(u,b) ] du

so the product Q*S is formed at each node and marginalised once. The
expectation of a product is not the product of the expectations: the two
differ by int_0^t Cov_b(Q,S) du, the association between quality of life
and survival that the joint model estimates.

The survival function is resolved from the fitted family with
merlin_p_getpf(), which returns a pointer to the prediction function for
whichever family the survival submodel was fitted with. Every survival
family merlin can predict from is therefore available, and merlin raises
its own error for any it cannot.
*/

// survival function of the survival submodel, resolved from its family.
// getpf dispatches on JM.model, so that has to be set before the call.
`PS' jmqaly_sfunc(`JM' JM)
{
	JM.model = JM.modtoind = strtoreal(st_global("jmqaly_stmodel"))
	return(merlin_p_getpf(JM,"survival"))
}

// quality of life on the natural (utility) scale
//
// merlin_util_expval() applies the submodel's inverse link, so this is
// correct for any family merlin supports an expected value for, not just
// gaussian.
//
// transform() covers the separate case of a continuous utility modelled on
// a transformed scale (e.g. gaussian on logit-transformed utilities): the
// link is the identity, so the back-transformation is applied here. It is
// applied to the mean trajectory, so what comes back is g^-1 of the linear
// predictor.
//
// The transform is passed in already resolved rather than read from the
// globals here: this is called once per quadrature node, and predictnl calls
// the whole user function again for every parameter it perturbs, so a string
// lookup in this function is repeated thousands of times per -jmqaly- call.
`RM' jmqaly_q(`JM' JM, `RC' t, `RS' dotrans, `RS' a, `RS' b)
{
	`RM' q

	q = merlin_util_expval(JM,t)
	if (!dotrans) return(q)
	return(a :* invlogit(q) :+ b)
}

// restricted QALYs: int_0^t Q(u) S(u) exp(-drate*u) du
`RM' jmqaly_qaly(`JM' JM, | `RC' t)
{
	`RS' ngq, lmod, smod, q, drate, ddisc, hasd
	`RS' dotrans, fmin, fmax, tmin, tmax, a, b
	`RM' gq, nodes, weights, yfit, sfit, res, dfac
	`RC' zi
	`PS' ps

	if (args()==1) t = merlin_util_timevar(JM)
	ngq 	= strtoreal(st_global("jmqaly_qintpoints"))
	gq 	= merlin_gq(ngq,"legendre")
	nodes 	= t :* J(rows(t),1,gq[,1]'):/2 :+ t:/2
	//weights on the natural scale: the survival is formed separately, so
	//there is nothing to be gained by carrying the weight in logs, and
	//log(0) at t=0 would be missing.
	weights = t :* J(rows(t),1,gq[,2]'):/2

	smod	= strtoreal(st_global("jmqaly_stmodel"))
	lmod 	= strtoreal(st_global("jmqaly_xtmodel"))
	ps	= jmqaly_sfunc(JM)

	//discount factor: exp(-r u) by default, (1+r)^-u with dmethod(discrete)
	drate	= strtoreal(st_global("jmqaly_drate"))
	if (missing(drate)) drate = 0
	ddisc	= st_global("jmqaly_ddisc")!=""
	hasd	= drate!=0

	//resolve the transform once, not once per quadrature node
	dotrans	= st_global("jmqaly_transform")!=""
	a = 1
	b = 0
	if (dotrans) {
		fmin = strtoreal(st_global("jmqaly_frommin"))
		fmax = strtoreal(st_global("jmqaly_frommax"))
		tmin = strtoreal(st_global("jmqaly_tomin"))
		tmax = strtoreal(st_global("jmqaly_tomax"))
		a    = (tmax - tmin) / (fmax - fmin)
		b    = tmax - a * fmax
	}

	res = 0
	for (q=1;q<=ngq;q++) {
		JM.model = JM.modtoind = lmod
		yfit = jmqaly_q(JM,nodes[,q],dotrans,a,b)
		JM.model = JM.modtoind = smod
		sfit = (*ps)(JM,nodes[,q])
		if (hasd) {
			if (ddisc)	dfac = (1+drate):^(-nodes[,q])
			else		dfac = exp(-drate :* nodes[,q])
			res = res :+ yfit :* sfit :* weights[,q] :* dfac
		}
		else	res = res :+ yfit :* sfit :* weights[,q]
	}
	//int_0^0 is zero whatever the integrand. merlin's rp basis is a spline
	//in log(t), so every quadrature node at t=0 is missing; the limit is
	//unambiguous here, so set it rather than propagating the missing.
	zi = selectindex(t:==0)
	if (rows(zi)) res[zi,] = J(rows(zi),cols(res),0)

	return(res)
}

// survival probability
`RM' jmqaly_surv(`JM' JM, | `RC' t)
{
	`PS' ps

	if (args()==1) t = merlin_util_timevar(JM)

	ps = jmqaly_sfunc(JM)
	JM.model = JM.modtoind = strtoreal(st_global("jmqaly_stmodel"))
	return((*ps)(JM,t))
}

// restricted mean survival time: int_0^t S(u) du
`RM' jmqaly_rmst(`JM' JM, | `RC' t)
{
	`RS' ngq, smod, q
	`RM' gq, nodes, weights, res
	`RC' zi
	`PS' ps

	if (args()==1) t = merlin_util_timevar(JM)
	ngq 	= strtoreal(st_global("jmqaly_qintpoints"))
	gq 	= merlin_gq(ngq,"legendre")
	nodes 	= t :* J(rows(t),1,gq[,1]'):/2 :+ t:/2
	weights = t :* J(rows(t),1,gq[,2]'):/2

	smod	= strtoreal(st_global("jmqaly_stmodel"))
	ps	= jmqaly_sfunc(JM)

	res = 0
	for (q=1;q<=ngq;q++) {
		JM.model = JM.modtoind = smod
		res = res :+ weights[,q] :* (*ps)(JM,nodes[,q])
	}
	//int_0^0 is zero whatever the integrand. merlin's rp basis is a spline
	//in log(t), so every quadrature node at t=0 is missing; the limit is
	//unambiguous here, so set it rather than propagating the missing.
	zi = selectindex(t:==0)
	if (rows(zi)) res[zi,] = J(rows(zi),cols(res),0)

	return(res)
}

end
