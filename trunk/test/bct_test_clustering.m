if ~exist("subtest", "var") || ~subtest
	source bct_test_setup.m
end

% clustering_coef_bd
for i = 1:size(m)(2)
	bct_test(sprintf("clustering_coef_bd %s", mname{i}), clustering_coef_bd(m{i}) == clustering_coef_bd_cpp(m{i})')
end

% clustering_coef_bu
for i = 1:size(m)(2)
	bct_test(sprintf("clustering_coef_bu %s", mname{i}), clustering_coef_bu(m{i}) == clustering_coef_bu_cpp(m{i})')
end

% clustering_coef_wd
for i = 1:size(m)(2)
	bct_test(sprintf("clustering_coef_wd %s", mname{i}), clustering_coef_wd(m{i}) == clustering_coef_wd_cpp(m{i})')
end

% clustering_coef_wu
for i = 1:size(m)(2)
	bct_test(sprintf("clustering_coef_wu %s", mname{i}), abs(clustering_coef_wu(m{i}) - clustering_coef_wu_cpp(m{i})') < 1e-6)
end

if ~exist("subtest", "var") || ~subtest
	printf("Failures: %d\n", failures)
	clear;
end
