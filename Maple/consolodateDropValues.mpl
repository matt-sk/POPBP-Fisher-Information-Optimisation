LIBRARY_FOLDER := "../lib";
read( cat(LIBRARY_FOLDER, "/lambda.mpl" ) ):

DATA_FOLDER := cat("Optimisation/n=",N,"/data");
DropValuesFile := cat( DATA_FOLDER, "/DropValues.m" );

for lambda in lambdaValues[N] do
	read( cat("dropValues (n=",N,",lambda=",lambda,").m") );
	d[lambda] := DV;
end do;

save( d, DropValuesFile );