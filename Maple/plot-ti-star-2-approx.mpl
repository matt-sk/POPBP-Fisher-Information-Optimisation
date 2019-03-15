LIBRARY_FOLDER := FileTools[JoinPath](["..", "lib"] ):
DATA_FOLDER := FileTools[JoinPath](["Optimisation", "n=2", "data"] ):

read FileTools[JoinPath]([LIBRARY_FOLDER, "opt.mpl"]):		# Optimisation routines (properly accounting for boundaries). Defines `opt`
read FileTools[JoinPath]([LIBRARY_FOLDER, "Fisher.mpl"]):	# External Fisher information routines. Defines `Fisher:-FI`
read FileTools[JoinPath]([LIBRARY_FOLDER, "lambda.mpl"]):	# Lambda values used for optimisation. Defines `lambda`

# Indexing function for tables which return NULL for unassigned entries.
`index/null` := proc(Idx::list,Tbl::table,Entry::list)
	if (nargs = 2) then
		if assigned(Tbl[op(Idx)]) then Tbl[op(Idx)];
		else NULL;
		end if;
	else
		Tbl[op(Idx)] := op(Entry);
	end if;
end proc: 


# Read the value of lambda to be used
N := 2: # FI_Approx is specific to the case of N=2.
lambda := lambdaValues[N][L];

# Initialise Timestamp
ID_String := sprintf( "n=%a, lambda=%a, timestamp=%s", N, lambda, StringTools[FormatTime]("%Y-%m-%d %X") );

optData := table():

# Function to perform optimisation, and save results into the external table.
local O := proc( p )
	local calc;
	global optData, opt;
	option remember;

	calc := opt[2]( p, lambda, Fisher:-FI_Approx[2] ):
	optData[p] := table( null, [ `FI*`=calc[1],`t1*`= calc[2][1], `t2*`=1.] ); # The `t2*` entry is entirely unnecessary, but kept for clarity.

	return calc;
end proc:


# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
# -= Find the drop value for this approximation  =-
# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
read "findDrop.mpl": # Utility to find the drop point in a plot of ti* (for Fisher Information).
infolevel['Fisher'] := 5;
DV := findDrop(N, lambda, Fisher:-FI_Approx[2], threshold = 6); 
minDV, maxDV := lhs(DV[1]), rhs(DV[1]);
midDV := (minDV+maxDV)/2;


# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
# -= Plot the curves for t1* as p varies         =-
# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-

# Part 1 is a straight line at `t1*`=1.
pl[1] := Matrix( [[0., 1.],[minDV,1.]] );
# We have not used O(), so optData is not updated automatically. We do so manually.
optData[0.] := table( null, [`FI*`=Fisher:-FI_Approx[2](1.,1.,0.,lambda), `t1*`=1., `t2*`=1.] ):
optData[minDV] := table( null, [`FI*`=Fisher:-FI_Approx[2](1.,1.,minDV,lambda), `t1*`=1., `t2*`=1.]):
pValues[`pt1`] := {0., minDV}:

# Plot FI* for the first part.
P := plot( p->Fisher:-FI_Approx[2](1.,1.,p,lambda), 0. .. minDV ):

pValues[`pt2`] := {}:

# For part 2 we plot the curve, and extract the matrix.
P := plot( p->O( p )[2][1], rhs(DV[1])..HFloat(1.) ):  # HFloat is needed for optData[1.] to correctly evaluate for pt[2] later. (Not sure why).
pl[2] := plottools[getdata](P)[3]; 

# The separator is vertical (at p=midDV) from t1*=1 down to whaterver t1* is at the first point of part 2 (row 1 column 2).
sep[1] := Matrix( [[midDV,1.],[midDV,pl[2][1][2]]] ):

plotData[`t1*`][`data`] := convert( pl, list );
plotData[`t1*`][`separators`] := convert( sep, list );


# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
# -= Construct the CSV file for pgfPlots         =-
# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
header := "p,fi*,t1*,s1":

# Function to take a table entry in the OptData table, and produce a sequence of the FI* and t1* falues. Useful in the Matrix construction below
flattenOptData := T-> <T[`FI*`] | T[`t1*`] | 'NULL'> :

# Convert the matrix data to columns.
for i to N do
	P := LinearAlgebra[Column]( pl[i], 1 );
	numPts := LinearAlgebra[Dimension]( P ):
	# Create a matrix of the OptData list. Each row is p, FI* and t1* in that order.
	pt[i] := Matrix( [seq([P[k],flattenOptData(optData[P[k]])], k=1..numPts)] ):
	CSVpt[i] := Export( pt[i], target=direct, format="CSV" ):
od:

# Hand code the separator (it's easier). We dont' need the fancy 4-line separator for N=2, but we'll keep it for consistency with the other files.
# Matrix: <<minDV|'NULL'|'NULL'|'NULL'|'NULL'>,<midDV|'NULL'|'NULL'|'NULL'|sep[1][1][2]>,<midDV|'NULL'|'NULL'|'NULL'|sep[1][2][2]>,<maxDV|'NULL'|'NULL'|'NULL'|'NULL'>>:
CSVsep[1] := sprintf( "%a,,,\n%a,,,%a\n%a,,,%a\n%a,,,\n", minDV, midDV, sep[1][1][2], midDV, sep[1][2][2], maxDV ):


# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
# -= Save results to files                       =-
# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-

# plotData matrices (for use with Maple)
mapleDataFileName := FileTools[JoinPath]([DATA_FOLDER, sprintf("ti*-approx plotData lambda=%.1f.m", lambda)]):
save( plotData, mapleDataFileName );

# CSV file (for use with pgfplots)
CSVFileName := FileTools[JoinPath]([DATA_FOLDER, sprintf("ti*-approx plotData lambda=%.1f.csv", lambda)]):
CSVfile := fopen( CSVFileName, WRITE ):
fprintf( CSVfile, "%s\n%s\n%s\n%s", header, CSVpt[1], CSVsep[1], CSVpt[2] ):
fclose( CSVfile ):
