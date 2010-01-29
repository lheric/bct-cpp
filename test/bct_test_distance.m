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

% distance_bin
for i = 1:size(m)(2)
	D=distance_bin(m{i});
	D_cpp = distance_bin_cpp(m{i});
	bct_test(sprintf("distance_bin %s", mname{i}), D == D_cpp);
end

% distance_wei
for i = 1:size(m)(2)
	D = distance_wei(m{i});
	D_cpp = distance_wei_cpp(m{i});
	bct_test(sprintf("distance_wei %s", mname{i}), D == D_cpp);
	
	%charpath ------
	[lambda,ecc,radius,diameter] = charpath(D);
	lambda_cpp = charpath_lambda_cpp(D);
	bct_test(sprintf("charpath lambda %s", mname{i}), lambda == lambda_cpp);
	[ecc_cpp, radius_cpp, diameter_cpp] = charpath_ecc_cpp(D);
	bct_test(sprintf("charpath ecc %s", mname{i}), ecc' == ecc_cpp);
	bct_test(sprintf("charpath radius %s", mname{i}), radius == radius_cpp);
	bct_test(sprintf("charpath diameter %s", mname{i}), diameter == diameter_cpp);
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

% findpaths and cycprob
sources = unique(floor(30*rand(1,5)+1));   % 5 random source nodes in the range (1, 30)
path_len_max = floor(rand()*2+2);  %random number in the range [2,4]
for i = 1:size(m)(2)
	[Pq,tpath,plq,qstop,allpths,util] = findpaths(m{i}, sources, path_len_max, 1);
	[allpaths_cpp Pq_cpp qstop_cpp] = findpaths_cpp(m{i}, sources, path_len_max, 1);
	bct_test(sprintf("findpaths allpaths %s", mname{i}), allpths == allpaths_cpp);
	bct_test(sprintf("findpaths Pq %s", mname{i}), Pq == Pq_cpp);
	bct_test(sprintf("findpaths qstop %s", mname{i}), qstop == qstop_cpp);
	
	%cycprob ---------
	[fcyc,pcyc] = cycprob(Pq);
	fcyc_cpp = cycprob_fcyc_cpp(Pq_cpp, qstop_cpp);
	pcyc_cpp = cycprob_pcyc_cpp(Pq_cpp, qstop_cpp);
	bct_test(sprintf("cycprob fcyc %s", mname{i}), fcyc == fcyc_cpp);
	bct_test(sprintf("cycprob pcyc %s", mname{i}), pcyc == pcyc_cpp);
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
	
	%charpath ------
	[lambda,ecc,radius,diameter] = charpath(D);
	lambda_cpp = charpath_lambda_cpp(D);
	bct_test(sprintf("charpath lambda %s", mname{i}), lambda == lambda_cpp);
	[ecc_cpp, radius_cpp, diameter_cpp] = charpath_ecc_cpp(D);
	bct_test(sprintf("charpath ecc %s", mname{i}), ecc' == ecc_cpp);
	bct_test(sprintf("charpath radius %s", mname{i}), radius == radius_cpp);
	bct_test(sprintf("charpath diameter %s", mname{i}), diameter == diameter_cpp);	
end

if ~exist("subtest", "var") || ~subtest
	printf("Failures: %d\n", failures)
	clear;
end
