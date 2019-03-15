# Set variables for important folders
LIBRARY_FOLDER := FileTools[JoinPath](["..", "lib"], base = worksheetdir):

DATA_FOLDER := FileTools[JoinPath](["Optimisation", cat("n=",N), "data"], base = worksheetdir):

# Read the necessary libraries.
read FileTools[JoinPath]([LIBRARY_FOLDER, "opt.mpl"]):		# Optimisation routines (properly accounting for boundaries). Defines `opt`
read FileTools[JoinPath]([LIBRARY_FOLDER, "Fisher.mpl"]):	# External Fisher information routines. Defines `Fisher:-FI`
read FileTools[JoinPath]([LIBRARY_FOLDER, "lambda.mpl"]):	# Lambda values used for optimisation. Defines `lambda`

# Read LAMBDA from passed in L value.
LAMBDA := lambdaValues[N][L];

# Initialise file names (this makes the code to use the files be a bit more readabile)
DropValuesFile := FileTools[JoinPath]( [DATA_FOLDER, "DropValues.m"] );
DataFileName := FileTools[JoinPath]( [DATA_FOLDER,sprintf("plotData lambda=%a.m", lambda)] );
SnapshotFileName := sprintf(".snapshot,n=%a,lambda=%a", N, LAMBDA );

# Initialise Timestamp
ID_String := sprintf( "n=%a, lambda=%a, timestamp=%s", N, LAMBDA, StringTools[FormatTime]("%Y-%m-%d %X") );

# Initialise values used in the computation.
delta := 0.005;	# The minimum distance between successive sampled values of p.

# Set up the calculation functions.
try
	read( SnapshotFileName );
catch :
	F := proc( P, L, N ) 
		option remember;
		local St,Fin, Opt;
		global Fisher;

		userinfo( 5, Fisher, `p=`||P );

		St:=time[real]();
		Opt:=opt[N]( P, L, Fisher:-FI[N] )[2];
		Fin:=time[real]();

		userinfo( 5, Fisher, (Fin-St), `seconds` );
		return Opt;
	end proc;
end try;

Snapshot_F := proc( )
	global F;
	local Opt;

	# Call F with the passed in paramters, and then save a snapshot of F, so we don't lose the computed (remembered) data before the plot is completed.
	Opt := F( _passed );
	save( F, SnapshotFileName );

	return Opt;
end proc;

# Set the drop values (including an extra drop value of 1.0..1.0 to make boundary checking easier.
read( DropValuesFile );
DV := [ op(d[LAMBDA]), 1.0..1.0 ];

# Sanity Check. Make sure that all dropvalues are less than 1, in which case we know that 
if false in map( i->is(rhs(i)<1.0), DV[1..-2] ) then
	ERROR( "Found a drop value not less than p=1.0" );
	quit;
end if:

# Known results for p=1. (This case performs poorly for the above code)
add((t[i]-t[i-1])^2/(exp(-lambda*t[i-1])-exp(-lambda*t[i])), i = 1 .. N); 
eval(%, {t[0] = 0, t[N] = 1, lambda=LAMBDA}); 
FI[p=1] := unapply(%, t, lambda);
# The `ti*`[p=1] function doesn't handle an input of 0 very well, so we make sure we optimise over a float very close to 0, up to 1.
# Also, the optimisation is specific to the value of N, so we need a piecewise to use the correct one.
# Finally, we only need to check the interior, because we know p=1 is after all the drop-values, so t_i^* has droped from 1 for all i. 
`ti*`[p=1] := piecewise(
	N=2,Optimization[NLPSolve]( t1->FI[p=1]([t1]), (1e-10)..1, maximize, method=branchandbound)[2],
	N=3,Optimization[NLPSolve]( (t1,t2)->FI[p=1]([t1,t2]), { (t1,t2)->t1-t2 }, {}, (1e-10)..1, (1e-10)..1, maximize, initialpoint=[1-1/sqrt(2), 1/sqrt(2)] )[2],
	N=4,Optimization[NLPSolve]( (t1,t2,t3)->FI[p=1]([t1,t2,t3]), { (t1,t2,t3)->t1-t2, (t1,t2,t3)->t2-t3 }, {}, (1e-10)..1, (1e-10)..1, (1e-10)..1, maximize, initialpoint=[1/4, 1-1/sqrt(2), 1/sqrt(2)] )[2]
);

# Hard-code the correct value for p=1.0 into the function F (because the c++ implementation doesn't handle this case)
F( 1.0, LAMBDA, Fisher:-FI[N] ) := `ti*`[p=1];

# Set the info level so we get diagnostic output suitable for recording times.
infolevel['Fisher']:=5;

# PrintLevel so that we may see the output of nested for loops.
printlevel := 3;

# Plot the parts of the graph, independently.
for i from 1 to N-1 do
	# We are producing the parts of the graph for ti*

	# Horozontal line before the drop point D_i(lambda)
	pl[i-1] := plot( <<0.0|1.0>,<lhs(DV[i])|1.0>> );

	# Plot the parts after the drop point D_i(lambda)
	for part from i to N-1 do
		StPt := rhs(DV[part]);
		EndPt := lhs(DV[part+1]);
		NumPts := max( ceil((EndPt-StPt)/delta), 3 );

		pl[part] := plot( P->Snapshot_F( P, LAMBDA, N )[i], StPt..EndPt, numpoints=NumPts );
	end do;
	# Extract the plot data matrices from the stored plots. (The data matrix is alwasy the 3rd element of the list returned by plottools[getdata]())
	pl := map( x->plottools[getdata](x)[3], pl );

	# Plot the vertical separation lines.
	for part from i to N-1 do;
		dv := (lhs(DV[part]) + rhs(DV[part]))/2;
		y := pl[part-1][-1][2], pl[part][1][2]; # NOTE: part numbers increment, so pl[0] becomes plotData[ti*][data][1], etc.
		sep[part] := plot( <<dv|y[1]>,<dv|y[2]>> );
	end do:
	# Extract the plot data matrices from the stored plots.
	sep := map( x->plottools[getdata](x)[3], sep );

	# Convert the plot data and separators data into lists.
	plotData[`t`||i||`*`][`data`] := convert( pl, list );
	plotData[`t`||i||`*`][`separators`] := convert( sep, list );

	# Unassign the tables, so that we don't save parts from previous iterations (each new i has fewer parts, starting from larger indeces).
	unassign( 'pl', 'sep' );
end do;

# Save the data, and remove the snapshot file.
save( plotData, DataFileName );
fremove( SnapshotFileName );
