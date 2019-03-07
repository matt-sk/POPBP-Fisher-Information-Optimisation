read( "findDrop.mpl" ):

LIBRARY_FOLDER := "../lib";
read( cat(LIBRARY_FOLDER, "/opt.mpl") ):	# Optimisation routines (properly accounting for boundaries). Defines `opt`
read( cat(LIBRARY_FOLDER, "/Fisher.mpl") ):	# External Fisher information routines. Defines `Fisher:-FI`
read( cat(LIBRARY_FOLDER, "/lambda.mpl") ):	# Lambda values used for optimisation. Defines `lambda`


storeDrop := proc( N::posint, lambda::realcons, FI::procedure )
	local ST, T;
	global DV, THRESHOLD;

	ST := time[real](): 
	if assigned( DV ) then 
		DV := findDrop(N, lambda, FI, bounds=DV, threshold = THRESHOLD): 
	else
		DV := findDrop(N, lambda, FI, threshold = THRESHOLD): 
	fi;
	T:=time[real]()-ST:

	print(lambda, DV, T):
end proc:


# The first time findDrop is called it takes a little longer than expected.
# this is likely due to some sort of internal setup. We run it once so that this doesn't affect our timings.
findDrop( N, 0.01, Fisher:-FI[N], threshold = 1 ):

# Read the value of lambda to be used
lambda := lambdaValues[N][L]:

# Set the output file (one file per combination of N and lambda)
outputFile := cat("dropValues (n=",N,",lambda=",lambda,").m" );

# Read the outputfile. This allows us to refine a previously computed range.
# If the read fails for any reason we quietly ignore it and proceed from scratch.
try 
	read( outputFile );
	if assigned('DV') then printf( "Continuing on from previous calculation: %a", op(DV) ); fi:
catch: 
end try;

# Find the drop value, save it in the variable `DV`, and save it to the output file.
infolevel['Fisher']:=5;
storeDrop( N, lambda, Fisher:-FI[N] );
save( DV, outputFile ):
