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
	bct_test(sprintf("betweenness_bin %s", mname{i}), all(betweenness_bin(m{i}) == betweenness_bin_cpp(m{i})))
end

% betweenness_wei
for i = 1:size(m)(2)
	bct_test(sprintf("betweenness_wei %s", mname{i}), all(betweenness_wei(m{i}) == betweenness_wei_cpp(m{i})'))
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

% findpaths
sources = unique(floor(30*rand(1,5)));
path_len_max = floor(rand()*2+2); %random number in the range [2,4]
for i = 1:size(m)(2)
	[Pq,tpath,plq,qstop,allpths,util] = findpaths(m{i}, sources, path_len_max, 1);
	allpaths_cpp = findpaths_cpp(m{i}, sources, path_len_max, 1);
	bct_test(sprintf("findpaths %s", mname{i}), all(allpths == allpaths_cpp));
	clear allpaths_cpp;
	clear Pq tpath plq qstop allpths util;
end

% jdegree
for i = 1:size(m)(2)
	[J J_od J_id J_bl] = jdegree(m{i});
	bct_test(sprintf("jdegree %s", mname{i}), all(J == jdegree_cpp(m{i})))
	bct_test(sprintf("jdegree_id %s", mname{i}), (J_id == jdegree_id_cpp(J)))
	bct_test(sprintf("jdegree_od %s", mname{i}), (J_od == jdegree_od_cpp(J)))
	bct_test(sprintf("jdegree_bl %s", mname{i}), (J_bl == jdegree_bl_cpp(J)))
end

% matching_ind
for i = 1:size(m)(2)
	[Min Mout Mall] = matching_ind(m{i});
	bct_test(sprintf("matching_ind %s", mname{i}), all(Mall == matching_ind_cpp(m{i})))
	%bct_test(sprintf("matching_ind_in %s", mname{i}), all(Min == matching_ind_in_cpp(m{i})))
	%bct_test(sprintf("matching_ind_out %s", mname{i}), all(Mout == matching_ind_out_cpp(m{i})))
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
