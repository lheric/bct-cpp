if ~exist("subtest", "var") || ~subtest
	source bct_test_setup.m
end

% breadth
for i = 1:size(m)(2)
	source = floor(length(m{i}) * rand()) + 1;
	[distance branch] = breadth(m{i}, source);
	[distance_cpp branch_cpp] = breadth_cpp(m{i}, source);
	bct_test(sprintf("breadth %s distance", mname{i}), distance == distance_cpp);
	bct_test(sprintf("breadth %s branch", mname{i}), branch == branch_cpp);
end

% breadthdist
for i = 1:size(m)(2)
	[R D] = breadthdist(m{i});
	[R_cpp D_cpp] = breadthdist_cpp(m{i});
	bct_test(sprintf("breadthdist %s R", mname{i}), R == R_cpp);
	bct_test(sprintf("breadthdist %s D", mname{i}), D == D_cpp);
end

% charpath
for i = 1:size(m)(2)
	D = breadthdist_cpp(m{i});
	[lambda ecc radius diameter] = charpath(D);
	[ecc_cpp radius_cpp diameter_cpp] = charpath_ecc_cpp(D);
	bct_test(sprintf("charpath %s lambda", mname{i}), lambda == charpath_lambda_cpp(D));
	bct_test(sprintf("charpath %s ecc", mname{i}), ecc == ecc_cpp');
	bct_test(sprintf("charpath %s radius", mname{i}), radius == radius_cpp);
	bct_test(sprintf("charpath %s diameter", mname{i}), diameter == diameter_cpp);
end

% cycprob
for i = 1:size(m)(2)
	sources = unique(floor(length(m{i}) * rand(1, 5))) + 1;
	qmax = 3;
	Pq = findpaths_cpp(m{i}, sources, qmax);
	[fcyc pcyc] = cycprob(Pq);
	bct_test(sprintf("cycprob %s fcyc", mname{i}), fcyc == cycprob_fcyc_cpp(Pq));
	bct_test(sprintf("cycprob %s pcyc", mname{i}), pcyc == cycprob_pcyc_cpp(Pq));
end

% distance_bin
for i = 1:size(m)(2)
	bct_test(sprintf("distance_bin %s", mname{i}), distance_bin(m{i}) == distance_bin_cpp(m{i}));
end

% distance_wei
for i = 1:size(m)(2)
	bct_test(sprintf("distance_wei %s", mname{i}), distance_wei(m{i}) == distance_wei_cpp(m{i}));
end

% efficiency
for i = 1:size(m)(2)
	E_local = efficiency(m{i},1);
	E_local_cpp = efficiency_local_cpp(m{i});
	E_global = efficiency(m{i});
	E_global_cpp = efficiency_global_cpp(m{i});
	bct_test(sprintf("efficiency_local %s", mname{i}), abs(E_local - E_local_cpp) < 1e-6);
	bct_test(sprintf("efficiency_global %s", mname{i}), abs(E_global - E_global_cpp) < 1e-6);
end

% findpaths
for i = 1:size(m)(2)
	sources = unique(floor(length(m{i}) * rand(1, 5))) + 1;
	qmax = 3;
	[Pq tpath plq qstop allpths util] = findpaths(m{i}, sources, qmax, 1);
	[Pq_cpp tpath_cpp plq_cpp qstop_cpp allpths_cpp util_cpp] = findpaths_cpp(m{i}, sources, qmax);
	bct_test(sprintf("findpaths %s Pq", mname{i}), Pq == Pq_cpp);
	bct_test(sprintf("findpaths %s tpath", mname{i}), tpath == tpath_cpp);
	bct_test(sprintf("findpaths %s plq", mname{i}), plq == plq_cpp);
	bct_test(sprintf("findpaths %s qstop", mname{i}), qstop == qstop_cpp);
	bct_test(sprintf("findpaths %s allpths", mname{i}), allpths == allpths_cpp);
	bct_test(sprintf("findpaths %s util", mname{i}), util == util_cpp);
end

% findwalks
for i = 1:size(m)(2)
	[Wq,twalk,wlq] = findwalks(m{i});
	N = size(m{i},1);
	[wlq_cpp,twalk_cpp,Wq_cpp]  = findwalks_cpp(m{i}, N); 
	%the second parameter is needed for testing from matlab but may not be necessary when used elsewhere
	bct_test(sprintf("findwalks wlq %s", mname{i}), wlq == wlq_cpp);
	bct_test(sprintf("findwalks Wq %s", mname{i}), Wq == Wq_cpp);
	bct_test(sprintf("findwalks twalk %s", mname{i}), twalk == twalk_cpp);
end

%reachdist
for i = 1:size(m)(2)
	[R,D] = reachdist(m{i});
	D_cpp = reachdist_cpp(m{i});
	bct_test(sprintf("reachdist %s", mname{i}), D == D_cpp);
end

if ~exist("subtest", "var") || ~subtest
	printf("Failures: %d\n", failures)
	clear;
end
