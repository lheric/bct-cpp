load cat
cat_all = CIJall;
cat_ctx = CIJctx;

load fve30
fve30 = CIJ;

load fve32
fve32 = CIJ;

load macaque47
macaque47 = CIJ;

load macaque71
macaque71 = CIJ;

m = {cat_all cat_ctx fve30 fve32 macaque47 macaque71};
mname = {"cat_all" "cat_ctx" "fve30" "fve32" "macaque47" "macaque71"};

global failures
failures = 0;

% assortativity
for i = 1:size(m)(2)
	bct_test(sprintf("assortativity_und %s", mname{i}), assortativity(m{i}, 0) == assortativity_und_cpp(m{i}))
	bct_test(sprintf("assortativity_dir %s", mname{i}), assortativity(m{i}, 1) == assortativity_dir_cpp(m{i}))
end

% betweenness_bin
for i = 1:size(m)(2)
	[EBC BC] = edge_betweenness_bin(m{i});
	bct_test(sprintf("betweenness_bin %s", mname{i}), all(BC == betweenness_bin_cpp(m{i})'))
end

% betweenness_wei
for i = 1:size(m)(2)
	[EBC BC] = edge_betweenness_wei(m{i});
	bct_test(sprintf("betweenness_wei %s", mname{i}), all(BC == betweenness_wei_cpp(m{i})'))
end

% breadth
for i = 1:size(m)(2)
	source = floor(5*rand())+1;  %a random source node in the range (1,5)
	[distance, branch] = breadth(m{i}, source);
	distance_cpp = breadth_cpp(m{i}, source);
	bct_test(sprintf("breadth %s", mname{i}), all(distance == distance_cpp));
end

% breadthdist
for i = 1:size(m)(2)
	[R, D] = breadthdist(m{i});
	distance_cpp = breadthdist_cpp(m{i});
	bct_test(sprintf("breadthdist %s", mname{i}), all(D == distance_cpp));
end

% clustering_coef_bd
for i = 1:size(m)(2)
	bct_test(sprintf("clustering_coef_bd %s", mname{i}), all(clustering_coef_bd(m{i}) == clustering_coef_bd_cpp(m{i})'))
end

% clustering_coef_bu
for i = 1:size(m)(2)
	bct_test(sprintf("clustering_coef_bu %s", mname{i}), all(clustering_coef_bu(m{i}) == clustering_coef_bu_cpp(m{i})'))
end

% clustering_coef_wd
for i = 1:size(m)(2)
	bct_test(sprintf("clustering_coef_wd %s", mname{i}), all(clustering_coef_wd(m{i}) == clustering_coef_wd_cpp(m{i})'))
end

% clustering_coef_wu
for i = 1:size(m)(2)
	bct_test(sprintf("clustering_coef_wu %s", mname{i}), all(abs(clustering_coef_wu(m{i}) - clustering_coef_wu_cpp(m{i})') < 1e-6))
end

% degrees_dir
for i = 1:size(m)(2)
	[id od deg] = degrees_dir(m{i});
	bct_test(sprintf("degrees_dir %s", mname{i}), all(deg == degrees_dir_cpp(m{i})))
end

% degrees_und
for i = 1:size(m)(2)
	bct_test(sprintf("degrees_und %s", mname{i}), all(degrees_und(m{i}) == degrees_und_cpp(m{i})))
end

% density_dir
for i = 1:size(m)(2)
	bct_test(sprintf("density_dir %s", mname{i}), density_dir(m{i}) == density_dir_cpp(m{i}))
end

% density_und
for i = 1:size(m)(2)
	bct_test(sprintf("density_und %s", mname{i}), density_und(m{i}) == density_und_cpp(m{i}))
end

% distance_bin
for i = 1:size(m)(2)
	D=distance_bin(m{i});
	D_cpp = distance_bin_cpp(m{i});
	bct_test(sprintf("distance_bin %s", mname{i}), all(D == D_cpp));
end

% distance_wei
for i = 1:size(m)(2)
	D = distance_wei(m{i});
	D_cpp = distance_wei_cpp(m{i});
	bct_test(sprintf("distance_wei %s", mname{i}), all(D == D_cpp));
	
	%charpath ------
	[lambda,ecc,radius,diameter] = charpath(D);
	lambda_cpp = charpath_lambda_cpp(D);
	bct_test(sprintf("charpath lambda %s", mname{i}), lambda == lambda_cpp);
	[ecc_cpp, radius_cpp, diameter_cpp] = charpath_ecc_cpp(D);
	bct_test(sprintf("charpath ecc %s", mname{i}), all(ecc' == ecc_cpp));
	bct_test(sprintf("charpath radius %s", mname{i}), radius == radius_cpp);
	bct_test(sprintf("charpath diameter %s", mname{i}), diameter == diameter_cpp);
end

% efficiency
for i = 1:size(m)(2)
	E_local = efficiency(m{i},1);
	E_local_cpp = efficiency_local_cpp(m{i});
	E_global = efficiency(m{i});
	E_global_cpp = efficiency_global_cpp(m{i});
	bct_test(sprintf("efficiency_local %s", mname{i}), all(E_local - E_local_cpp)<0.0000001);
	bct_test(sprintf("efficiency_global %s", mname{i}), all(E_global - E_global_cpp)<0.0000001);
end

% edge_betweenness_bin
for i = 1:size(m)(2)
	bct_test(sprintf("edge_betweenness_bin %s", mname{i}), all(all(edge_betweenness_bin(m{i}) == edge_betweenness_bin_cpp(m{i}))))
end

% edge_betweenness_wei
for i = 1:size(m)(2)
	bct_test(sprintf("edge_betweenness_wei %s", mname{i}), all(all(edge_betweenness_wei(m{i}) == edge_betweenness_wei_cpp(m{i}))))
end

% erange
for i = 1:size(m)(2)
	bct_test(sprintf("erange %s", mname{i}), all(all(erange(m{i}) == erange_cpp(m{i}))))
end

% findpaths and cycprob
sources = unique(floor(30*rand(1,5)+1));   % 5 random source nodes in the range (1, 30)
path_len_max = floor(rand()*2+2);  %random number in the range [2,4]
disp(['sources' sources]); 
disp(['path_len_max' path_len_max]); 
for i = 1:size(m)(2)
	[Pq,tpath,plq,qstop,allpths,util] = findpaths(m{i}, sources, path_len_max, 1);
	[allpaths_cpp Pq_cpp qstop_cpp] = findpaths_cpp(m{i}, sources, path_len_max, 1);
	bct_test(sprintf("findpaths allpaths %s", mname{i}), all(allpths == allpaths_cpp));
	bct_test(sprintf("findpaths Pq %s", mname{i}), all(Pq == Pq_cpp));
	bct_test(sprintf("findpaths qstop %s", mname{i}), (qstop == qstop_cpp));
	
	%cycprob ---------
	[fcyc,pcyc] = cycprob(Pq);
	fcyc_cpp = cycprob_fcyc_cpp(Pq_cpp, qstop_cpp);
	pcyc_cpp = cycprob_pcyc_cpp(Pq_cpp, qstop_cpp);
	bct_test(sprintf("cycprob fcyc %s", mname{i}), all(fcyc == fcyc_cpp));
	bct_test(sprintf("cycprob pcyc %s", mname{i}), all(pcyc == pcyc_cpp));
end

% findwalks
for i = 1:size(m)(2)
	[Wq,twalk,wlq] = findwalks(m{i});
	N = size(m{i},1);
	[wlq_cpp,twalk_cpp,Wq_cpp]  = findwalks_cpp(m{i}, N); 
	%the second parameter is needed for testing from matlab but may not be necessary when used elsewhere
	bct_test(sprintf("findwalks wlq %s", mname{i}), all(wlq == wlq_cpp));
	bct_test(sprintf("findwalks Wq %s", mname{i}), all(Wq == Wq_cpp));
	bct_test(sprintf("findwalks twalk %s", mname{i}), twalk == twalk_cpp);
end

% jdegree
for i = 1:size(m)(2)
	[J J_od J_id J_bl] = jdegree(m{i});
	bct_test(sprintf("jdegree %s", mname{i}), all(J == jdegree_cpp(m{i})))
	bct_test(sprintf("jdegree_id %s", mname{i}), (J_id == jdegree_id_cpp(J)))
	bct_test(sprintf("jdegree_od %s", mname{i}), (J_od == jdegree_od_cpp(J)))
	bct_test(sprintf("jdegree_bl %s", mname{i}), (J_bl == jdegree_bl_cpp(J)))
end

% latmio_dir
for i = 1:size(m)(2)
	R = latmio_dir(m{i},2);
	deg_latmio_dir = degrees_dir(R);
	R_cpp = latmio_dir_cpp(m{i}, 2);
	deg_latmio_dir_cpp = degrees_dir(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("latmio_und_connected %s", mname{i}), all(abs(deg_latmio_dir-deg_latmio_dir_cpp) == 0)); 
end

% latmio_dir_connected
for i = 1:size(m)(2)
	R = latmio_dir_connected(m{i},2);
	deg_latmio_dir = degrees_dir(R);
	R_cpp = latmio_dir_connected_cpp(m{i}, 2);
	deg_latmio_dir_cpp = degrees_dir(R_cpp);
	bct_test(sprintf("latmio_dir_connected %s", mname{i}), all(abs(deg_latmio_dir-deg_latmio_dir_cpp) == 0));
end

% latmio_und_connected
for i = 1:size(m)(2)
	R = latmio_und_connected(m{i},2);
	deg_latmio_und = degrees_und(R);
	R_cpp = latmio_und_connected_cpp(m{i}, 2);
	deg_latmio_und_cpp = degrees_und(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("latmio_und_connected %s", mname{i}), all(abs(deg_latmio_und-deg_latmio_und_cpp) < 10)); 
end

% latmio_und
for i = 1:size(m)(2)
	R = latmio_und(m{i},2);
	deg_latmio_und = degrees_und(R);
	R_cpp = latmio_und_cpp(m{i}, 2);
	deg_latmio_und_cpp = degrees_und(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("latmio_und_connected %s", mname{i}), all(abs(deg_latmio_und-deg_latmio_und_cpp) < 10)); 
end


% matching_ind
for i = 1:size(m)(2)
	[Min Mout Mall] = matching_ind(m{i});
	bct_test(sprintf("matching_ind %s", mname{i}), all(Mall == matching_ind_cpp(m{i})))
end	

%reachdist
for i = 1:size(m)(2)
	[R,D] = reachdist(m{i});
	D_cpp = reachdist_cpp(m{i});
	bct_test(sprintf("reachdist %s", mname{i}), all(D == D_cpp));
	
	%charpath ------
	[lambda,ecc,radius,diameter] = charpath(D);
	lambda_cpp = charpath_lambda_cpp(D);
	bct_test(sprintf("charpath lambda %s", mname{i}), lambda == lambda_cpp);
	[ecc_cpp, radius_cpp, diameter_cpp] = charpath_ecc_cpp(D);
	bct_test(sprintf("charpath ecc %s", mname{i}), all(ecc' == ecc_cpp));
	bct_test(sprintf("charpath radius %s", mname{i}), radius == radius_cpp);
	bct_test(sprintf("charpath diameter %s", mname{i}), diameter == diameter_cpp);	
end

% strengths_und
for i = 1:size(m)(2)
	bct_test(sprintf("strengths_und %s", mname{i}), all(strengths_und(m{i}) == strengths_und_cpp(m{i})))
end

% strengths_dir
for i = 1:size(m)(2)
	[is os str] = strengths_dir(m{i});
	bct_test(sprintf("strengths_dir %s", mname{i}), all(str == strengths_dir_cpp(m{i})))
end

printf("Failures: %d\n", failures)

clear;
