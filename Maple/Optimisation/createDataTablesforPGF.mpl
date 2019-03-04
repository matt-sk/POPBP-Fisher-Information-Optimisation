# We extract the data for the ti* plots, from the matrices produced in their computation, and save it in text format.

# Ideally, this would have been done as part of the creation of the plots, wherein we could also have made sure that 
# all curves for a given part used the same p values, but neither idea occured before computing the graphs and the time 
# requried for computation makes re-computation prohibitive. As such, we extract and manipulate the data here with some 
# small inconvenience.

with( LinearAlgebra ):

LIBRARY_FOLDER := FileTools[JoinPath](["..", "..", "lib"] ):

read FileTools[JoinPath]([LIBRARY_FOLDER, "lambda.mpl"]):

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
	# Initialise the data folder for this N.
	DATA_FOLDER := FileTools[JoinPath]( [cat("n=",N), "data"] ):

	# We make sure that the variables p, tk*, and sk are all unset so that the header (below) can be constructed.
	unassign( 'p' ); #, seq(cat(t, k, `*`), k = 1 .. N-1), seq(cat(s, k), k = 1 .. N-1) );

	# Construct the header that will be used for all graph data for this value of N.
	header := cat( "p", seq(cat(",t",i,"*"), i=1..N-1), seq(cat(",s",i), i=1..N-1) ):

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

		# We break the data up into parts, breaking up lines that span parts into separate sub-lines traversing the same parts.

		# We divide the range of p into parts and transitions. Parts are the values of p between the bounds of subsequent drop values.
		# Transitions are the values of p within the bounds of drop points. 
		# Parts contain only data plots, and Transitiosn might contain data or separator plots.

		# -= =-=-=-=-=-=-=-=-=-=-=-=-= =-
		# -= Special case for part 1.  =-
		# -= =-=-=-=-=-=-=-=-=-=-=-=-= =-
		# All curves for `tj*` have the same data (`tj*`=1). The plot data for this part is stored in the first matrix of t1*
		p := LinearAlgebra[Column]( plotData[`t1*`][data][1], 1 );
		`t*` := LinearAlgebra[Column]( plotData[`t1*`][data][1], 2 );
		s := <'NULL', 'NULL'>;

		# Note that the plot data for part 1 for curves t2*, ... t(n-1)* is stored in an odd way.
		# Since the curves for those variables continue as a striaght line through the separator, 
		# and into the next part(s), the entire straight line is stored as a single 'part' in 
		# P[2],...,P[N-1] respectively that spans several parts and separator spaces.
		# We separate this data out into each individual part and transition,
		# The data for part 1 is all the same for both the ti* curves, as well as the separtor curves.
		partData[1] := Matrix( [ p, seq( `t*`, k=1..(N-1) ), seq( s, k=1..(N-1) ) ] );

		# Clean up.
		unassign( 'p', '`t*`', 's' );

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
				s[j]	:= piecewise( j > tr,  <'NULL','NULL','NULL','NULL'>, <'NULL', LinearAlgebra[Column]( plotData[t||j||`*`][separators][-(N-tr)],2), 'NULL'>);
			end do;
			
			# The order of columns in the matrix is `p`, `t1*`, ... `tN*`, `s1`, ... , `sN`
			transitionData[tr] := Matrix( [ p, seq( `t*`[k], k=1..(N-1) ), seq( s[k], k=1..(N-1) ) ] );
		end do;
		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=- =-

		# Cleanup Vectors that will be reused.
		unassign( 'p', '`t*`', 's' );

		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-= =-
		# -= Process parts (i.e., values of p within the bounds of  \(\mathfrak{D}_i\)) =-
		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-= =-
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

			# Now we can construct the matrix for part i
			numPts := nops(pValues);
			p := convert( sort(pValues), Vector );
			s := Vector( numPts, k->NULL );
			for j to N-1 do
				`t*`[j]	:= Vector( numPts, k->data[j][p[k]] );
			end do;

			partData[pt] := Matrix( [ p, seq( `t*`[k], k=1..(N-1) ), seq( s, k=1..(N-1) ) ] );
		end do;
		# -= =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-= =-

		# Cleanup Vectors that will be reused.
		unassign( 'p', '`t*`', 's' );

		# -= =-=-=-=-=-=-=-=-=-=-= =-
		# -= Create the CSV files
		# -= =-=-=-=-=-=-=-=-=-=-= =-


		# Create the format string for sprintf. Each partData[k] and transitionData[k] is a separate 
		# block of the data file Separated by an empty line.
		# We always have the header, N parts and N-1 separators, interleaved. In total that is N+1 strings.
		formatStr := cat("%s\n","%s"$N);
		for k to N-1 do
		   strData[k] := sprintf( "%s\n%s\n", Export( partData[k], target=direct, format="CSV"), Export( transitionData[k], target=direct, format="CSV"));
		end do;
		strData[N] := Export( partData[k], target=direct, format="CSV");

		CSVfileName := FileTools[JoinPath]([DATA_FOLDER, sprintf("plotData lambda=%.1f.csv", lambda)]);
		CSVfile := fopen( CSVfileName, WRITE);
		dataCSV := fprintf( CSVfile, formatStr, header, seq(strData[k], k=1..N) );
		fclose( CSVfile );
		# -= =-=-=-=-=-=-=-=-=-=-= =-
	end do 
end do: