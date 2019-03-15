try
	# Check for needed assignments.
	if not assigned( LIBRARY_FOLDER ) then
		LIBRARY_FOLDER := ".":
	end if:

	if not assigned( NUMTHREADS ) then
		ERROR( "NUMTHREADS not assigned" ):
	end if:

	Fisher := module( )
		export FI, FI_Approx;
		local external, external_threaded, _cpp_library;

		# Initialiation
		FI := table();
		external := table();
		external_threaded := table();

		_cpp_library := cat( LIBRARY_FOLDER, "/Fisher.so" );

		# Define the external functions, both threaded and unthreaded.
		external[2] := define_external( "FI2", t1::(float[8]), t2::(float[8]), p::(float[8]), lambda::(float[8]), RETURN::(float[8]), LIB = _cpp_library ):
		external[3] := define_external( "FI3", t1::(float[8]), t2::(float[8]), t3::(float[8]), p::(float[8]), lambda::(float[8]), RETURN::(float[8]), LIB = _cpp_library ):
		external[4] := define_external( "FI4", t1::(float[8]), t2::(float[8]), t3::(float[8]), t4::(float[8]), p::(float[8]), lambda::(float[8]), RETURN::(float[8]), LIB = _cpp_library ):
		external[5] := define_external( "FI5", t1::(float[8]), t2::(float[8]), t3::(float[8]), t4::(float[8]), t5::(float[8]), p::(float[8]), lambda::(float[8]), RETURN::(float[8]), LIB = _cpp_library ):

		external_threaded[2] := define_external( "FI2_threaded", t1::(float[8]), t2::(float[8]), p::(float[8]), lambda::(float[8]), _numThreads::(integer[4]), RETURN::(float[8]), LIB = _cpp_library ): 
		external_threaded[3] := define_external( "FI3_threaded", t1::(float[8]), t2::(float[8]), t3::(float[8]), p::(float[8]), lambda::(float[8]), _numThreads::(integer[4]), RETURN::(float[8]), LIB = _cpp_library ): 
		external_threaded[4] := define_external( "FI4_threaded", t1::(float[8]), t2::(float[8]), t3::(float[8]), t4::(float[8]), p::(float[8]), lambda::(float[8]), _numThreads::(integer[4]), RETURN::(float[8]), LIB = _cpp_library ): 
		external_threaded[5] := define_external( "FI5_threaded", t1::(float[8]), t2::(float[8]), t3::(float[8]), t4::(float[8]), t5::(float[8]), p::(float[8]), lambda::(float[8]), _numThreads::(integer[4]), RETURN::(float[8]), LIB = _cpp_library ): 

		# Define the functions for calculating Fisher Information.
		# Note that FI[2] is not threaded, but the rest of the external.
		FI[2] := (t1, t2, p, lambda) -> external[2](t1, t2, p, lambda ): 
		FI[3] := (t1, t2, t3, p, lambda) -> external_threaded[3](t1, t2, t3, p, lambda, NUMTHREADS ): 
		FI[4] := (t1, t2, t3, t4, p, lambda) -> external_threaded[4](t1, t2, t3, t4, p, lambda, NUMTHREADS ): 
		FI[5] := (t1, t2, t3, t4, t5, p, lambda) -> external_threaded[5](t1, t2, t3, t4, t5, p, lambda, NUMTHREADS ): 

		# The approximation to FI[2]. 
		FI_Approx[2] := (1+p/v[1])*p*(p+(1-p)*(p*v[1, 2]+(1-p)*v[2])-(1-p)*(p*v[1, 2]+(1-p)*v[2])^2)*((t2-t1)*p+(1-p)*t2*v[1])^2/((p+(1-p)*v[1])^2*(p+p*(1-p)*v[1, 2]+(1-p)^2*v[2])^2*(1-p*v[1, 2]-(1-p)*v[2]))-p/(p+(1-p)*v[1])*(p*(t2-t1)^2*(p+(1-p)*(1-v[1, 2])*v[1, 2])/((1-v[1, 2])*(p+(1-p)*v[1, 2])^2))+p*t1^2*(p+(1-p)*(1-v[1])*v[1])/((1-v[1])*(p+(1-p)*v[1])^2):
		FI_Approx[2] := eval(FI_Approx[2], {v[1] = exp(-lambda*t1), v[2] = exp(-lambda*t2), v[1, 2] = exp(-lambda*(t2-t1))}):
		# Note that there are undefined values when t1=0 and/or t1=t2, so we must take the limits and use a piecewise function to choose between the different expressions.
		FI_Approx[2] := piecewise(
			t1=0.,			piecewise( t2=0., limit( limit( FI_Approx[2], t1=0 ), t2=0 ), limit( FI_Approx[2], t1=0 ) ),
			t2-t1=0.,		limit( FI_Approx[2], t2=t1 ), # t1 and t2 are the same, up to numeric precision.
			FI_Approx[2] # Default
		):
		FI_Approx[2] := unapply( FI_Approx[2], t1, t2, p, lambda ):
	end module;

catch:
	# Mimick the output of an error call before quitting with an error value of 1.
	str := "":
	if lastexception[1] <> 0 then str := sprintf( " (in %a)", lastexception[1] ) fi:
	printf( "Error,%s %s\n", str, StringTools[FormatMessage](lastexception[2..-1]) ):

end try: