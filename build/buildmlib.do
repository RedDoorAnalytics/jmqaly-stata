// build ljmqaly.mlib from ./mata
//
// run from the repository root:
//	do build/buildmlib.do
//
// The library is replaced by `mlib create, replace' at the END, once every
// source below has compiled, and is never erased first: that way a compile
// error leaves the previous library standing rather than none at all, so the
// next command to run fails on the compile error rather than on a missing
// function.

mata: mata set matastrict on
qui {
	do "./mata/jmqaly.mata"

	mata: mata mlib create ljmqaly, dir(.) replace
	mata: mata mlib add    ljmqaly jmqaly_*(), dir(.)
	mata: mata d *()
	mata mata clear
	mata mata mlib index

	// the library is compiled strict; restore the default so a user's own
	// mata code is not forced to satisfy matastrict
	mata: mata set matastrict off
}

di as result "ljmqaly.mlib built"
