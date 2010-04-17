bct_test_setup

load motif34lib

% motif3generate
[M ID N] = motif3generate_cpp;
bct_test("motif3generate M3", M3 == M)
bct_test("motif3generate ID3", ID3 == ID')
bct_test("motif3generate N3", N3 == N')

% motif3funct_bin
for i = 1:size(m)(2)
	[f F] = motif3funct_bin(m{i});
	[f_cpp F_cpp] = motif3funct_bin_cpp(m{i});
	bct_test(sprintf("motif3funct_bin %s f", mname{i}), f == f_cpp')
	bct_test(sprintf("motif3funct_bin %s F", mname{i}), F == F_cpp)
end

% motif3funct_wei
for i = 1:size(m)(2)
	[I Q F] = motif3funct_wei(m{i});
	[I_cpp Q_cpp F_cpp] = motif3funct_wei_cpp(m{i});
	bct_test(sprintf("motif3funct_wei %s I", mname{i}), I == I_cpp)
	bct_test(sprintf("motif3funct_wei %s Q", mname{i}), Q == Q_cpp)
	bct_test(sprintf("motif3funct_wei %s F", mname{i}), F == F_cpp)
end

% motif3struct_bin
for i = 1:size(m)(2)
	[f F] = motif3struct_bin(m{i});
	[f_cpp F_cpp] = motif3struct_bin_cpp(m{i});
	bct_test(sprintf("motif3struct_bin %s f", mname{i}), f == f_cpp')
	bct_test(sprintf("motif3struct_bin %s F", mname{i}), F == F_cpp)
end

% motif3struct_wei
for i = 1:size(m)(2)
	[I Q F] = motif3struct_wei(m{i});
	[I_cpp Q_cpp F_cpp] = motif3struct_wei_cpp(m{i});
	bct_test(sprintf("motif3struct_wei %s I", mname{i}), I == I_cpp)
	bct_test(sprintf("motif3struct_wei %s Q", mname{i}), Q == Q_cpp)
	bct_test(sprintf("motif3struct_wei %s F", mname{i}), F == F_cpp)
end

% motif4generate
[M ID N] = motif4generate_cpp;
bct_test("motif4generate M4", M4 == M)
bct_test("motif4generate ID4", ID4 == ID')
bct_test("motif4generate N4", N4 == N')

% motif4funct_bin
for i = 3:size(m)(2)
	[f F] = motif4funct_bin(m{i});
	[f_cpp F_cpp] = motif4funct_bin_cpp(m{i});
	bct_test(sprintf("motif4funct_bin %s f", mname{i}), f == f_cpp')
	bct_test(sprintf("motif4funct_bin %s F", mname{i}), F == F_cpp)
end

% motif4struct_bin
for i = 1:size(m)(2)
	[f F] = motif4struct_bin(m{i});
	[f_cpp F_cpp] = motif4struct_bin_cpp(m{i});
	bct_test(sprintf("motif4struct_bin %s f", mname{i}), f == f_cpp')
	bct_test(sprintf("motif4struct_bin %s F", mname{i}), F == F_cpp)
end

% motif4struct_wei
for i = 1:size(m)(2)
	[I Q F] = motif4struct_wei(m{i});
	[I_cpp Q_cpp F_cpp] = motif4struct_wei_cpp(m{i});
	bct_test(sprintf("motif4struct_wei %s I", mname{i}), I == I_cpp)
	bct_test(sprintf("motif4struct_wei %s Q", mname{i}), Q == Q_cpp)
	bct_test(sprintf("motif4struct_wei %s F", mname{i}), F == F_cpp)
end

bct_test_teardown
