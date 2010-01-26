if ~exist("subtest", "var") || ~subtest
	source bct_test_setup.m
end

% betweenness_bin
for i = 1:size(m)(2)
	[EBC BC] = edge_betweenness_bin(m{i});
	bct_test(sprintf("betweenness_bin %s", mname{i}), BC == betweenness_bin_cpp(m{i})')
end

% betweenness_wei
for i = 1:size(m)(2)
	[EBC BC] = edge_betweenness_wei(m{i});
	bct_test(sprintf("betweenness_wei %s", mname{i}), BC == betweenness_wei_cpp(m{i})')
end

% edge_betweenness_bin
for i = 1:size(m)(2)
	bct_test(sprintf("edge_betweenness_bin %s", mname{i}), edge_betweenness_bin(m{i}) == edge_betweenness_bin_cpp(m{i}))
end

% edge_betweenness_wei
for i = 1:size(m)(2)
	bct_test(sprintf("edge_betweenness_wei %s", mname{i}), edge_betweenness_wei(m{i}) == edge_betweenness_wei_cpp(m{i}))
end

% erange
for i = 1:size(m)(2)
	bct_test(sprintf("erange %s", mname{i}), erange(m{i}) == erange_cpp(m{i}))
end

if ~exist("subtest", "var") || ~subtest
	printf("Failures: %d\n", failures)
	clear;
end
