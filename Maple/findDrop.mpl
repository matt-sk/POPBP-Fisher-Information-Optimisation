findDrop := proc (N::posint, lambda, FI::procedure, { bounds::~Array(range(float[8])) := [(HFloat(0)..HFloat(1))$(N-1)], threshold::posint := Digits-1 }) 
	option hfloat;
	local chk, testP, minP, maxP, P, i, j, ST, T;

	# Sanity Checking
	if ArrayDims(bounds) <> 1 .. (N-1) then 
		ERROR( "expected keyword parameter bounds to have array dimensions [1 .. %2], but received array with dimensions %3", procname, N-1, [ArrayDims(bounds)] );
	end if; 

	# Initialise
	minP := map( lhs, bounds ); #Array( 1..N, lhs(bounds) );
	maxP := map( rhs, bounds ); #Array( 1..N, rhs(bounds) );

	# We find the drop point for each t_i^* except for t_N^* which we know must always be 1, so has no drop point.
	for i to N-1 do
		# For each i we perform a binary search until the min and max values our threshold of each other.
		while 10^(-threshold) < abs(minP[i]-maxP[i]) do 
			# Choose the midpoint between min and max P values, and optimise for that new P.
			testP := (1/2)*minP[i]+(1/2)*maxP[i];
			ST := time[real](); chk := opt[N](testP, lambda, FI)[2]; T := time[real]()-ST;
			
			# Print out the time taken for record purposes.
			userinfo( 5, Fisher, `Optimisation time for N=` || N || ` lambda=` || lambda || ` p=` || testP || `: `, `` || T || ` seconds` );

			# Check to see if t_i^* is equal to 1 (up to precision) and adjust min or max accordingly.
			if abs(chk[i]-1.0) < 10^(-Digits+1) then minP[i] := testP else maxP[i] := testP end if;

			# Update min and max for subsequent i values based on what we've seen so far. 
			# (This should speed up future iterations of the i-loop)
			for j from i+1 to N-1 do 
				if 10^(-Digits+1) <= abs(chk[j]-1.0) then 
					maxP[j] := min(maxP[j], testP) 
				else 
					minP[j] := max(minP[j], testP) 
				end if 
			end do;
		end do;

		# We know that the drop value is between min and max, and that both are within the required threshold of each other. 
		P[i] := minP[i]..maxP[i];
	end do;

	# Return the ranges as a list.
	return convert( P, list );
end proc;
