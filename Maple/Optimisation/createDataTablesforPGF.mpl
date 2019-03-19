# We extract the data for the ti* plots, from the matrices produced in their computation, and save it in text format.

# Ideally, this would have been done as part of the creation of the plots, wherein we could also have made sure that 
# all curves for a given part used the same p values, but neither idea occured before computing the graphs and the time 
# requried for computation makes re-computation prohibitive. As such, we extract and manipulate the data here with some 
# small inconvenience.

with( LinearAlgebra ):

LIBRARY_FOLDER := FileTools[JoinPath](["..", "..", "lib"] ):

read FileTools[JoinPath]([LIBRARY_FOLDER, "Fisher.mpl"]):	# External Fisher information routines. Defines `Fisher:-FI`
read FileTools[JoinPath]([LIBRARY_FOLDER, "lambda.mpl"]):	# Lambda values used for optimisation. Defines `lambdaValues`

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

# Initialise Timestamp
ID_String := sprintf( "createDataTablesforPGF, timestamp=%s", StringTools[FormatTime]("%Y-%m-%d %X") );


# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
# -= Main Loop: extract and format the data  =-
# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
for N from 2 to 4 do 
	printf( "n=%a:", N ):
	# Initialise the data folder for this N.
	DATA_FOLDER := FileTools[JoinPath]( [cat("n=",N), "data"] ):

	# Known results for p=1 based on SBP. (This case performs poorly for the C++ code)
	unassign( 't', 'lambda' ):
	add((t[i]-t[i-1])^2/(exp(-lambda*t[i-1])-exp(-lambda*t[i])), i = 1 .. N); 
	eval(%, {t[0] = 0, t[N] = 1}); 
	FI[SBP] := unapply(%, t, lambda);

	# Construct the header that will be used for all graph data for this value of N.
	header := cat( "p", ",FI*", seq(cat(",t",i,"*"), i=1..N-1), seq(cat(",s",i), i=1..N-1) ):

	# Timing information.
	START_N := time[real]():

	# Produce the graphs for each value of lambda.
	for lambda in lambdaValues[N] do 
		try
			# Read the data file for n and lambda.
			mapleDataFileName := FileTools[JoinPath]([DATA_FOLDER, sprintf("plotData lambda=%a.m", lambda)]);
			read mapleDataFileName; 
		catch : 
			# Output error and continue to next possibility.
			WARNING( "couldn't read file for n=%1, lambda=%2", N, lambda );
			next; 
		end try;

		# Timing information.
		START_lambda := time[real]():

		# We break the data up into parts, breaking up curves that span multiple parts into separate sub-curves each 
		# traversing one of the those parts.

		# We divide the range of p into parts and transitions. Parts are the values of p between the bounds of subsequent drop values.
		# Transitions are the values of p within the bounds of drop points. 
		# Parts contain only data plots, and Transitiosn might contain data or separator plots.

		# -= =-=-=-=-=-=-=-=-=-=-=-=-= =-
		# -= Special case for part 1.  =-
		# -= =-=-=-=-=-=-=-=-=-=-=-=-= =-
		# All curves for `ti*` have the same data (`ti*`=1). The plot data for this part is stored in the first matrix of t1*
		p := Column( plotData[`t1*`][data][1], 1 );

		# Plot FI* Data for this part.
		PlotData	:= plot( p->Fisher:-FI[N]( (1.0$N), p, lambda ), p[1]..p[2] );
		PlotData	:= plottools[getdata]( PlotData )[3]:

		# Extract data from the plot, and construct the relevant entires for `ti*` and `s` variables.
		p,`FI*`		:= Column( PlotData, 1 ), Column( PlotData, 2 ):

		numPts		:= Dimension( p ):
		s			:= Matrix( numPts, N-1, 'NULL' ): 
		`t*`		:= Matrix( numPts, N-1, (k,i)->piecewise(k=1 or k=numPts, 1, 'NULL') ): 

		# Note that the plot data for part 1 for curves t2*, ... t(n-1)* is stored in an odd way.
		# Since the curves for those variables continue as a striaght line through the separator, 
		# and into the next part(s), the entire straight line is stored as a single 'part' in 
		# P[2],...,P[N-1] respectively that spans several parts and separator spaces.
		# We separate this data out into each individual part and transition,
		# The data for part 1 is all the same for both the ti* curves, as well as the separtor curves.
		partData[1] := Matrix( [ p, `FI*`, `t*`, s ] );

		# Clean up.
		unassign( 'PlotData', 'p', '`FI*`', '`t*`', 's' );

		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=- =-
		# -= Process transitions (i.e., values of p within the bounds of  \(\mathfrak{D}_i\))  =-
		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=- =-
		for tr to N-1 do
			# Processing transition number tr

			# The min value of \(\mathfrak{D}_i\) is the x value of the last point in the first data matrix of t(tr)*.
			minDV := plotData[t||tr||`*`][data][1][-1][1];

			# The max value of \(\mathfrak{D}_i\) is the x value of the first point in the second data matrix of t(tr)*
			maxDV := plotData[t||tr||`*`][data][2][1][1];

			# The mid value of \(\mathfrak{D}_i\) is the x value of any point in the first separator matrix of t(tr)*
			midDV := plotData[t||tr||`*`][separators][1][1][1];

			# The p values will be minDVi, midDVi, midDVi, maxDVi. (Note the doubling of midDVi to allow a vertical separator line).
			p := <minDV, midDV, midDV, maxDV>;

			# For j > tr, `tj*` has data in transition tr, and `sj` has no separator in transition tr.
			# For j ≤ tr  `tj*` has no data in transition tr, and `sj` has a separator in transition tr.
			for j from 1 to N-1 do 
				`t*`[j]	:= piecewise( j > tr,  <1,'NULL','NULL',1>,  <'NULL','NULL','NULL','NULL'> );
				s[j]	:= piecewise( j > tr,  <'NULL','NULL','NULL','NULL'>, <'NULL', Column( plotData[t||j||`*`][separators][-(N-tr)],2), 'NULL'>);
			end do;

			# No FI* data is computed inside of a transition.
			`FI*` :=  <'NULL','NULL','NULL','NULL'>:
			
			# The order of columns in the matrix is `p`, `t1*`, ... `tN*`, `s1`, ... , `sN`
			transitionData[tr] := Matrix( [ p, `FI*`, seq( `t*`[k], k=1..(N-1) ), seq( s[k], k=1..(N-1) ) ] );
		end do;
		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=- =-

		# Cleanup Vectors that will be reused.
		unassign( 'p', '`FI*`', '`t*`', 's' );

		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
		# -= Process parts (i.e., values of p between the bounds of  \(\mathfrak{D}_i\)) =-
		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
		for pt from 2 to N do
			# Process part number pt

			# Initialise pValues as the empty set.
			pValues := {};

			for j to pt-1 do
				# Initillise data[j] as a table which defaults to NULL for unused indices.
				data[j] := table( null ); 

				# Plot data for tj* curve in Part i is stored as the [-(N-i)]th matrix for tj* part data.
				tempData := plotData[t||j||`*`][data][-(N-pt+1)];

				# Extract the individual p values for later use.
				pValues := pValues union convert( Column(tempData, 1), set );

				# Re-encode the tempData as a table with key p-value, and value tj*
				numPts := RowDimension( tempData );
				for k to numPts do
					# data[j][p] stores the value of tj* for probability p
					# Each row of tempData is <p | tj*> and we extract the data accordingly.
					data[j][tempData[k][1]] := tempData[k][2];
				end do;
			end do;

			# ti* for i > pt is still 1. We set these valeus now.
			for i from pt to N-1 do
				# Initillise data[j] as a table which defaults to NULL for unused indices.
				data[i] := table( null ); 

				data[i][min(pValues)] := 1.;
				data[i][max(pValues)] := 1.;
			end do;

			# Now we can construct the matrix for part pt
			numPts	:= nops(pValues);
			p		:= convert( sort(pValues), Vector );
			s		:= Matrix( numPts, N-1, 'NULL' );
			`t*`	:= Matrix( numPts, N-1, (k,i)->data[i][p[k]] ); # Columns correspond to `ti*`, rows to p-values.

			# Compute FI* for each p. (
			# Note: `t*[k]` is the kth row of the `t*` matrix contating t1*...t(pt-1)*
			`FI*` := Vector( numPts, 'NULL' ):
			for k to numPts do
				try:
					t := convert(`t*`[k][1..pt-1], list):
					`FI*`[k] := Fisher:-FI[N]( op(t), 1.0$(N-pt+1), p[k], lambda ):
				catch: # Do nothing if there's an error, just move on.
				end try:
			end do:

			# Concatenate the data for the part.
			partData[pt] := Matrix( [p, `FI*`, `t*`, s] ):
		end do;
		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-

		# Adjust for p=1. (Which should be the final row of the last part)
		`FI*`[-1] := FI[SBP]( convert( `t*`[-1], list ), lambda ):
		partData[N] := Matrix( [p, `FI*`, `t*`, s] ):

		# Timing information.
		END_lambda := time[real]():

		# Cleanup Vectors that will be reused.
		unassign( 'p', '`FI*`', '`t*`', 's' ):

		# -= =-=-=-=-=-=-=-=-=-=-= =-
		# -= Create the CSV files
		# -= =-=-=-=-=-=-=-=-=-=-= =-


		# Create the format string for sprintf. Each partData[k] and transitionData[k] is a separate 
		# block of the data file Separated by an empty line.
		# We always have the header, N parts and N-1 separators, interleaved. In total that is N+1 strings.
		formatStr := cat("%s\n","%s"$N):
		for k to N-1 do
		   strData[k] := sprintf( "%s\n%s\n", Export( partData[k], target=direct, format="CSV"), Export( transitionData[k], target=direct, format="CSV")):
		end do;
		strData[N] := Export( partData[k], target=direct, format="CSV"):

		CSVfileName := FileTools[JoinPath]([DATA_FOLDER, sprintf("plotData lambda=%.1f.csv", lambda)]):
		CSVfile := fopen( CSVfileName, WRITE):
		dataCSV := fprintf( CSVfile, formatStr, header, seq(strData[k], k=1..N) ):
		fclose( CSVfile ):
		# -= =-=-=-=-=-=-=-=-=-=-= =-

		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
		# -= Create and export the Maple plot matrices =-
		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= =-
		# Extract the first two columns (p and `FI*`) from the part matrices, and convert the resulting table into a list.
		plotData := convert( map( M->M[.., 1..2], partData), list ): # Index of “..” gives all entries.
		dataFileName := FileTools[JoinPath]([DATA_FOLDER, sprintf("FI* plotData lambda=%.1f.m", lambda)]):
		save( ID_String, plotData, dataFileName ):
		# -= =-=-=-=-=-=-=-=-=-=-= =-

		printf( " 𝜆=%.1g (%.1fs)", lambda, END_lambda - START_lambda ):
	end do:

	# Timing information.
	END_N := time[real]():

	printf( ". Total: %.1fs\n", END_N - START_N ):
end do: