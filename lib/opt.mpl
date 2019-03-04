opt[2] := proc (P, lambda, FI::procedure) 
	local maxVal, maxArg, tmp;
	uses Optimization;

	maxArg := [-infinity, -infinity];
	maxVal := -infinity;

	tmp := NLPSolve( (x) -> FI(x, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound);
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], 1.0] end if;

	return [maxVal, maxArg] 
end proc;

opt[3] := proc (P, lambda, FI::procedure) 
	local maxVal, maxArg, tmp;
	uses Optimization;

	maxArg := [-infinity, -infinity, -infinity];
	maxVal := -infinity;

	# Only if p > rhs(DV_2) ???
	tmp := NLPSolve( (x,y) -> FI(x, y, 1.0, P, lambda), { (x,y) -> x-y }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], 1.0] end if;

	tmp := NLPSolve( (x) -> FI(x, x, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], 1.0] end if;

	# Only if p >= lhs(DV_1) ???
	tmp := NLPSolve( (x) -> FI(x, 1.0, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], 1.0, 1.0] end if;
	# Not if p >= rhs(DV_1)
	tmp := FI( 1.0, 1.0, 1.0, P, lambda );
	if maxVal < tmp then maxVal, maxArg := tmp, [1.0, 1.0, 1.0] end if;

	return [maxVal, maxArg]
end proc;

opt[4] := proc (P, lambda, FI::procedure) 
	local maxVal, maxArg, tmp;
	uses Optimization;

	maxArg := [-infinity, -infinity, -infinity, -infinity];
	maxVal := -infinity;

	# Interior
	tmp := NLPSolve( (t1, t2, t3) -> FI(t1, t2, t3, 1.0, P, lambda), { (t1, t2, t3) -> t1-t2, (t1, t2, t3) -> t2-t3 }, {}, 0 .. 1, 0 .. 1, 0 .. 1, maximize, initialpoint = [.25, .25, .5] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], tmp[2][3], 1.0] end if;

	# 2-Dimeinsional Boundaries
	tmp := NLPSolve( (t1, t2) -> FI(t1, t2, 1.0, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], 1.0, 1.0] end if;

	tmp := NLPSolve( (t1, t2) -> FI(t1, t2, t2, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], tmp[2][2], 1.0] end if;

	tmp := NLPSolve( (t1, t2) -> FI(t1, t1, t2, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], tmp[2][2], 1.0] end if;

	# 1-Dimensional Boundaries
	tmp := NLPSolve( (t1) -> FI(t1, t1, t1, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], 1.0, 1.0, 1.0] end if;
	
	tmp := NLPSolve( (t1) -> FI(t1, t1, 1.0, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], 1.0, 1.0, 1.0] end if;
	
	tmp := NLPSolve( (t1) -> FI(t1, 1.0, 1.0, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], 1.0, 1.0, 1.0] end if;
	
	# 0-Dimensional Boundaries
	tmp := FI( 1.0, 1.0, 1.0, 1.0, P, lambda );
	if maxVal < tmp then maxVal, maxArg := tmp, [1.0, 1.0, 1.0, 1.0] end if;
	
	return [maxVal, maxArg] 
end proc;

opt[5] := proc (P, lambda, FI::procedure) 
	local maxVal, maxArg, tmp;
	uses Optimization;

	maxArg := [-infinity, -infinity, -infinity, -infinity, -infinity];
	maxVal := -infinity;

	# Interior
	tmp := NLPSolve( (t1, t2, t3, t4) -> FI(t1, t2, t3, t4, 1.0, P, lambda), { (t1, t2, t3, t4) -> t1-t2, (t1, t2, t3, t4) -> t2-t3, (t1, t2, t3, t4) -> t3-t4 }, {}, 0 .. 1, 0 .. 1, 0 .. 1, 0..1, maximize, initialpoint = [.25, .25, .25, .5] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], tmp[2][3], tmp[2][4], 1.0] end if;

	# 3-Dimensional Boundaries
	tmp := NLPSolve( (t1, t2, t3) -> FI(t1, t2, t3, 1.0, 1.0, P, lambda), { (t1, t2, t3) -> t1-t2, (t1, t2, t3) -> t2-t3 }, {}, 0 .. 1, 0 .. 1, 0 .. 1, maximize, initialpoint = [0.25, 0.25, 0.5] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], tmp[2][3], 1.0, 1.0] end if;

	tmp := NLPSolve( (t1, t2, t3) -> FI(t1, t1, t2, t3, 1.0, P, lambda), { (t1, t2, t3) -> t1-t2, (t1, t2, t3) -> t2-t3 }, {}, 0 .. 1, 0 .. 1, 0 .. 1, maximize, initialpoint = [0.25, 0.25, 0.5] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], tmp[2][2], tmp[2][3], 1.0] end if;

	tmp := NLPSolve( (t1, t2, t3) -> FI(t1, t2, t2, t3, 1.0, P, lambda), { (t1, t2, t3) -> t1-t2, (t1, t2, t3) -> t2-t3 }, {}, 0 .. 1, 0 .. 1, 0 .. 1, maximize, initialpoint = [0.25, 0.25, 0.5] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], tmp[2][2], tmp[2][3], 1.0] end if;

	tmp := NLPSolve( (t1, t2, t3) -> FI(t1, t2, t3, t3, 1.0, P, lambda), { (t1, t2, t3) -> t1-t2, (t1, t2, t3) -> t2-t3 }, {}, 0 .. 1, 0 .. 1, 0 .. 1, maximize, initialpoint = [0.25, 0.25, 0.5] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], tmp[2][3], tmp[2][3], 1.0] end if;

	# 2-Dimensional Boundaries
	tmp := NLPSolve( (t1, t2) -> FI(t1, t2, 1.0, 1.0, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], 1.0, 1.0, 1.0] end if;

	tmp := NLPSolve( (t1, t2) -> FI(t1, t1, t2, 1.0, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], tmp[2][2], 1.0, 1.0] end if;

	tmp := NLPSolve( (t1, t2) -> FI(t1, t2, t2, 1.0, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], tmp[2][2], 1.0, 1.0] end if;

	tmp := NLPSolve( (t1, t2) -> FI(t1, t1, t1, t2, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], tmp[2][1], tmp[2][2], 1.0] end if;

	tmp := NLPSolve( (t1, t2) -> FI(t1, t2, t2, t2, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][2], tmp[2][2], tmp[2][2], 1.0] end if;

	tmp := NLPSolve( (t1, t2) -> FI(t1, t1, t2, t2, 1.0, P, lambda), { (t1, t2) -> t1-t2 }, {}, 0 .. 1, 0 .. 1, maximize, initialpoint = [1-1/sqrt(2), 1/sqrt(2)] );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], tmp[2][2], tmp[2][2], 1.0] end if;

	# 1-Dimensional Boundaries
	tmp := NLPSolve( (t1) -> FI(t1, t1, t1, t1, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], tmp[2][1], tmp[2][1], 1.0] end if;
	
	tmp := NLPSolve( (t1) -> FI(t1, t1, t1, 1.0, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], tmp[2][1], 1.0, 1.0] end if;
	
	tmp := NLPSolve( (t1) -> FI(t1, t1, 1.0, 1.0, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], tmp[2][1], 1.0, 1.0, 1.0] end if;
	
	tmp := NLPSolve( (t1) -> FI(t1, 1.0, 1.0, 1.0, 1.0, P, lambda), {}, {}, 0 .. 1, maximize, method = branchandbound );
	if maxVal < tmp[1] then maxVal, maxArg := tmp[1], [tmp[2][1], 1.0, 1.0, 1.0, 1.0] end if;
	
	# 0-Dimensional Boundaries
	tmp := FI( 1.0, 1.0, 1.0, 1.0, 1.0, P, lambda );
	if maxVal < tmp then maxVal, maxArg := tmp, [1.0, 1.0, 1.0, 1.0, 1.0] end if;
	
	return [maxVal, maxArg] 
end proc;

