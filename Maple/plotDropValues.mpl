LIBRARY_FOLDER := FileTools[JoinPath](["..", "lib"] ):
DATA_FOLDER		:= FileTools[JoinPath](["Optimisation", cat("n=",N), "data"] ):

read FileTools[JoinPath]([LIBRARY_FOLDER, "opt.mpl"]):		# Optimisation routines (properly accounting for boundaries). Defines `opt`
read FileTools[JoinPath]([LIBRARY_FOLDER, "Fisher.mpl"]):	# External Fisher information routines. Defines `Fisher:-FI`

read "findDrop.mpl": # Utility function for finding the drop points (interval) for a given n, p and lambda. Defines `findDrop`

SnapshotFileName := cat( ".snapshot.plotDropValues n=",N,".m" );

PlotData := table():

# Set up the calculation functions.
F := proc( N, L ) 
	option remember;
	local St,Fin,DV;
	global Fisher, PlotData;

	userinfo( 5, Fisher, `p=`||P );

	St:=time[real]();
	DV:=findDrop( N, L, Fisher:-FI[N], threshold=3 );
	Fin:=time[real]();

	PlotData[L] := DV;
	save( PlotData, SnapshotFileName );

	userinfo( 4, Fisher, `lambda=`, L, `DV=`, DV, (Fin-St), `seconds` );
	return DV;
end proc:

# See if we are resuming from a previous computation or not,
try
	# Read the snapshot file and populate the remember table of the function.
	read( SnapshotFileName );

	lprint( "Resuming from previous computation." );
	# Initialise the remember table for F from the snapshot file data.
	for lambda in indices(PlotData, nolist) do
		F(N,lambda) := PlotData[lambda];
	end do;
	unassign( 'lambda' );
catch :
	# We are not resuming (we are computing from scratch).
	# We initialise the remember table for F with the values for lambda=1 assuming that the drop values are exactly 1.
	F(N,0.0) := [(1..1)$(N-1)];
	PlotData[0.0] := F(N,0.0);
end try:

# Set infolevel to 4 so that the info from F is printed out, but the info from findDrop is not.
infolevel[Fisher] := 4:


# Plot the bounds of the drop values.
# Note that these plots are simply to casue the relevant data for a smooth plot to be saved in both
# the PlotData table, and the remember table for F. We regenerate these plots later with more data.
for i from 1 to N-1 do
	P[i] := plot( [L->lhs(F(N,L)[i]),L->rhs(F(N,L)[i])], 0.0..5.0 )
end do:

# Colate the calculated data into a text table for use with LaTeX pgfplots.
# Note that due to the nature of the adaptive algorithm, it might be the case that
# P[i] and P[j] might contain data for some different lambda values whem j≠i.
# Each individually will have enough data to produce a smooth plot.
# The entire data (for all used values of lambda) is stored in the PlotData table.
# It seems a shame not to use the full data in all the plots, so we collate it all now.
with(LinearAlgebra):

# Extract the indices from the saved PlotData table. Note that these will be lists of [N,L].
idx := sort( [indices(PlotData, nolist)] ):

# Extract the lambda data from the indices, storing it as a column Vector.
lambdaValues := Vector(idx):

COL := plots[setcolours]():
for i from 1 to N-1 do
	# Extract the LHS and RHS of each of the ranges for \( \mathfrak{D}_i \)
	LHData[i] := Vector( nops(idx), k->lhs(PlotData[idx[k]][i]) ):
	RHData[i] := Vector( nops(idx), k->rhs(PlotData[idx[k]][i]) ):

	# Use this data to re-generate the plot of DV_i
	pl[1] := plot( <lambdaValues|LHData[i]>, colour=COL[i] ):
	pl[2] := plot( <lambdaValues|RHData[i]>, colour=COL[i] ):

	# Extract the plot data matrices from the stored plots. (The data matrix is alwasy the 3rd element of the list returned by plottools[getdata]())
	pl := map( x->plottools[getdata](x)[3], pl );
	plotData[i] := convert(pl, list): # Note that plotData and PlotData are different variables. (Capital vs lower case p).

	# Concatenate the minDV and maxDV data into a single matrix using the column vectors as, well, columns. 
	# This if for conversion to CSV later on.
	Data[i] := <LHData[i] | RHData[i]>: 
end do:;


# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
# -= Construct the CSV file for pgfPlots         =-
# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-

# Header row.
header := cat( "lambda", seq(cat(",min D",i,",max D",i), i=1..N-1) ):

# Express the dropValue data as a matrix, and then convert it into CSV.
CSVdata := Export( Matrix( [lambdaValues,seq(Data[k], k=1..N-1)] ), target=direct, format="CSV" ):


# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
# -= Save results to files                       =-
# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-

# plotData matrices (for use with Maple)
DataFileName := FileTools[JoinPath]( [DATA_FOLDER,"dropValues plotData.m"] );
save( plotData, DataFileName ):

# CSV file (for use with pgfplots)
CSVFileName := FileTools[JoinPath]([DATA_FOLDER, "dropValues plotData.csv"]):
CSVFile := fopen( CSVFileName, WRITE ):
fprintf( CSVFile, "%s\n%s", header, CSVdata ):
fclose( CSVFile ):
