try
	# Check for needed assignments.
	if not assigned( LIBRARY_FOLDER ) then
		LIBRARY_FOLDER := ".":
	end if:

	if not assigned( NUMTHREADS ) then
		ERROR( "NUMTHREADS not assigned" ):
	end if:

	Fisher := module( )
		export FI;
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
	end module;

catch:
	# Mimick the output of an error call before quitting with an error value of 1.
	str := "":
	if lastexception[1] <> 0 then str := sprintf( " (in %a)", lastexception[1] ) fi:
	printf( "Error,%s %s\n", str, StringTools[FormatMessage](lastexception[2..-1]) ):

end try: