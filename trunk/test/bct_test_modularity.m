if ~exist("subtest", "var") || ~subtest
	source bct_test_setup.m
end

if ~exist("subtest", "var") || ~subtest
	printf("Failures: %d\n", failures)
	clear;
end
