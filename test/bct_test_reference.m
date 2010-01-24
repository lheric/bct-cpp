if ~exist("subtest", "var") || ~subtest
	source bct_test_setup.m
end

% latmio_dir
for i = 1:size(m)(2)
	R = latmio_dir(m{i},2);
	deg_latmio_dir = degrees_dir(R);
	R_cpp = latmio_dir_cpp(m{i}, 2);
	deg_latmio_dir_cpp = degrees_dir(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("latmio_dir %s", mname{i}), deg_latmio_dir == deg_latmio_dir_cpp);
end

% latmio_dir_connected
for i = 1:size(m)(2)
	R = latmio_dir_connected(m{i},2);
	deg_latmio_dir = degrees_dir(R);
	R_cpp = latmio_dir_connected_cpp(m{i}, 2);
	deg_latmio_dir_cpp = degrees_dir(R_cpp);
	bct_test(sprintf("latmio_dir_connected %s", mname{i}), deg_latmio_dir == deg_latmio_dir_cpp);
end

% latmio_und_connected
for i = 1:size(m)(2)
	R = latmio_und_connected(m{i},2);
	deg_latmio_und = degrees_und(R);
	R_cpp = latmio_und_connected_cpp(m{i}, 2);
	deg_latmio_und_cpp = degrees_und(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("latmio_und_connected %s", mname{i}), abs(deg_latmio_und - deg_latmio_und_cpp) < 10);
end

% latmio_und
for i = 1:size(m)(2)
	R = latmio_und(m{i},2);
	deg_latmio_und = degrees_und(R);
	R_cpp = latmio_und_cpp(m{i}, 2);
	deg_latmio_und_cpp = degrees_und(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("latmio_und %s", mname{i}), abs(deg_latmio_und - deg_latmio_und_cpp) < 10);
end

% randmio_dir_connected
for i = 1:size(m)(2)
	R = randmio_dir_connected(m{i},2);
	deg_randmio_dir = degrees_dir(R);
	R_cpp = randmio_dir_connected_cpp(m{i}, 2);
	deg_randmio_dir_cpp = degrees_dir(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("randmio_dir_connected %s", mname{i}), deg_randmio_dir == deg_randmio_dir_cpp);
end

% randmio_dir
for i = 1:size(m)(2)
	R = randmio_dir(m{i},2);
	deg_randmio_dir = degrees_und(R);
	R_cpp = randmio_dir_cpp(m{i}, 2);
	deg_randmio_dir_cpp = degrees_dir(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("randmio_dir %s", mname{i}), deg_randmio_dir == deg_randmio_dir_cpp);
end

% randmio_und_connected
for i = 1:size(m)(2)
	R = randmio_und_connected(m{i},2);
	deg_randmio_und = degrees_und(R);
	R_cpp = randmio_und_connected_cpp(m{i}, 2);
	deg_randmio_und_cpp = degrees_und(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("randmio_und_connected %s", mname{i}), abs(deg_randmio_und - deg_randmio_und_cpp) < 10);
end

% randmio_und
for i = 1:size(m)(2)
	R = randmio_und(m{i},2);
	deg_randmio_und = degrees_und(R);
	R_cpp = randmio_und_cpp(m{i}, 2);
	deg_randmio_und_cpp = degrees_und(R_cpp);
	%this is only an approximate test condition
	bct_test(sprintf("randmio_und %s", mname{i}), abs(deg_randmio_und - deg_randmio_und_cpp) < 10);
end

if ~exist("subtest", "var") || ~subtest
	printf("Failures: %d\n", failures)
	clear;
end
