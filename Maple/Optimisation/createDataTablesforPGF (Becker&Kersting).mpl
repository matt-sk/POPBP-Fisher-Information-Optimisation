restart;
Digits := 30;

with(Optimization):

# The Becker and Kersting formula to calculate t[i]^%H; for the case where p = 1; appears to be an approximation (asympotic in n), not an exact value.
# Definition: Let `t*`[i]; denote the optimal value of t[i]; as calculated by the Becker and Kersting formula.
# Definition: Let `s*`[i]; denote the optimal value of t[i]; as calculated by the NLPSolve(); function.
# We will demonstrate that the Fisher information computed using conjugate(s[i]); is (slightly) greater than the Fisher information computed using conjugate(t[i]);.

# Basic function setup.

# Formula for Fisher Information for SBP (or, equivalently, a POSBP with p=1).
# NOTE: We explicitly set t[n] = 1 in this formula knowing that `t*[n]` must always be 1. As such this isn't the general Fisher information for that n, 
#       but it is sufficient for calculating the optimal values.
FI := proc( n )
	local EXP;

	# Take the basic expression and subsititute in t[0] = 0 and t[n] = 1 
	EXP := add((t[i] - t[i - 1])^2/(exp(-lambda*t[i - 1]) - exp(-lambda*t[i])), i = 1 .. n);
	EXP := eval(EXP, {t[0] = 0, t[n] = 1});

	# Turn this expression into a function of t and lambda, and return it.
	unapply( EXP, t, lambda );
end proc:

# Formula for Becker and Kersting approximation to ti* for the above fisher information function.
BKApprox := n -> unapply( subs( tau=1, [ seq( 3/lambda * log(1 + (i/n)*(exp(lambda*tau/3)-1)), i=1..n-1 ) ] ), lambda ):

OPT[2] := (lambda, IP) -> NLPSolve(t1 -> FI(2)([t1], lambda), 0.01 .. 0.99, maximize, method = branchandbound, evaluationlimit = 75, nodelimit = 75):
OPT[3] := (lambda, IP) -> NLPSolve((t1, t2) -> FI(3)([t1, t2], lambda), {(t1, t2) -> t1 - t2 + 0.01}, {}, 0.01 .. 0.99, 0.01 .. 0.99, maximize, initialpoint = IP):
OPT[4] := (lambda, IP) -> NLPSolve((t1, t2, t3) -> FI(4)([t1, t2, t3], lambda), {(t1, t2, t3) -> t1 - t2 + 0.01, (t1, t2, t3) -> t2 - t3 + 0.01}, {}, 0.01 .. 0.99, 0.01 .. 0.99, 0.01 .. 0.99, maximize, initialpoint = IP):

COMPARE := n -> proc(lambda) 
	local `CALC_t*`, `CALC_FI`, `OPT_t*`, `OPT_FI`; 
	global `FI`, `BKApprox`, `FI_p1`, `t*`, OPT;
	option remember;

	# CALC_?? uses the Becker and Kersting formulae.
	`CALC_t*` := BKApprox(n)(lambda);
	`CALC_FI` := FI(n)(`CALC_t*`, lambda);

	# OPT_?? uses numeric optimisation. We start at the poitn calculated by Becker and Kersting.
	`OPT_FI`, `OPT_t*` := op( OPT[n](lambda, `CALC_t*`) ); # NLPSolve returns optimal value first, then list of optimal paramaters

	table([
		CALC = table( [`t*` = `CALC_t*`, `FI` = `CALC_FI`] ),
		OPT  = table( [`t*` = `OPT_t*`,  `FI` = `OPT_FI`]  )
	]);
end proc:

WRITE_CSV := proc( COMPARISON_PLOT, { n::posint := NULL } ) 
	local `FI Comparison Data` := plottools[getdata](COMPARISON_PLOT)[1   , 3];
	local `t* Comparison Data` := plottools[getdata](COMPARISON_PLOT)[2.. , 3];

	local BASE_LAMBDA := LinearAlgebra[Column](`FI Comparison Data`, 1):
	local COMPARISON_LAMBDA, DATA;

	local CSV_HEADER := cat( "𝜆", ",FI*", seq(cat(",t",i,"*"), i=1..n-1) );

	local DATA_FOLDER := FileTools[JoinPath]( [cat("n=",N), "data"] ):
	local CSVfileName := FileTools[JoinPath]( [DATA_FOLDER, "plotData Backer&Kersting.csv"] ):
	local CSVfile, dataCSV;

	# Sanity Check. Make sure the lambda values are the same for all data.
	for COMPARISON_LAMBDA in map( LinearAlgebra[Column], {`t* Comparison Data`}, 1) do
		if not LinearAlgebra[Equal](BASE_LAMBDA, COMPARISON_LAMBDA ) then ERROR( "Lambda values differ in between FI and ti*" ) fi;
	end do;

	# Read the `t*` data column (column 2) from the `t* Comparison Data` matrices, and append to the `FI Comparison Data` matrix.
	DATA := Matrix( [`FI Comparison Data`, op(map( LinearAlgebra[Column], {`t* Comparison Data`}, 2))] );

	# Convert the matrix to CSV, and add the header.
	DATA := sprintf( "%s\n%s", CSV_HEADER, Export(DATA, target=direct, format="CSV") );

	# Write the data
	CSVfile := fopen( CSVfileName, WRITE );
	dataCSV := fprintf( CSVfile, "%s", DATA );
	fclose( CSVfile ):

	return DATA;
end proc:

# The case of n = 2. 
N := 2;
compare := COMPARE(N);

`t1* Comparison` := lambda -> compare(lambda)[`OPT`][`t*`][1] - compare(lambda)[`CALC`][`t*`][1];

# The Fisher Information can get quite large, and so the difference between the two can get quite large. We'll look at a percentage difference of the two that will be positive so long as OPTIMISE is greater, and negative so long as CALCULATED is greater.
`FI Comparison` := lambda -> 1 - compare(lambda)[`CALC`][FI]/compare(lambda)[`OPT`][FI];

COMPARISON_PLOT := plot( [ `FI Comparison`, `t1* Comparison` ] , 0 .. 50, numpoints = 1000 );

WRITE_CSV( COMPARISON_PLOT, n = N );

# The case of n = 3.
N := 3;
compare := COMPARE(N);

`t1* Comparison` := lambda -> compare(lambda)[`OPT`][`t*`][1] - compare(lambda)[`CALC`][`t*`][1];
`t2* Comparison` := lambda -> compare(lambda)[`OPT`][`t*`][2] - compare(lambda)[`CALC`][`t*`][2];

# The Fisher Information can get quite large, and so the difference between the two can get quite large. We'll look at a percentage difference of the two that will be positive so long as OPTIMISE is greater, and negative so long as CALCULATED is greater.
`FI Comparison` := lambda -> 1 - compare(lambda)[`CALC`][FI]/compare(lambda)[`OPT`][FI];

COMPARISON_PLOT := plot( [ `FI Comparison`, `t1* Comparison`, `t2* Comparison` ] , 0 .. 50, numpoints = 1000 );

WRITE_CSV( COMPARISON_PLOT, n = N );

# The case of n = 4.
N := 4;
compare := COMPARE(N);

`t1* Comparison` := lambda -> compare(lambda)[`OPT`][`t*`][1] - compare(lambda)[`CALC`][`t*`][1];
`t2* Comparison` := lambda -> compare(lambda)[`OPT`][`t*`][2] - compare(lambda)[`CALC`][`t*`][2];
`t3* Comparison` := lambda -> compare(lambda)[`OPT`][`t*`][3] - compare(lambda)[`CALC`][`t*`][3];

# The Fisher Information can get quite large, and so the difference between the two can get quite large. We'll look at a percentage difference of the two that will be positive so long as OPTIMISE is greater, and negative so long as CALCULATED is greater.
`FI Comparison` := lambda -> 1 - compare(lambda)[`CALC`][FI]/compare(lambda)[`OPT`][FI];

COMPARISON_PLOT := plot( [ `FI Comparison`, `t1* Comparison`, `t2* Comparison`, `t3* Comparison` ] , 0 .. 50, numpoints = 1000 );

WRITE_CSV( COMPARISON_PLOT, n = N );
