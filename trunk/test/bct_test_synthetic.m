bct_test_setup

N = floor(91 * rand()) + 10;
Kdir = floor((N * N - N + 1) * rand());
Kund = floor(Kdir / 2);

% makerandCIJ_bd
CIJ = makerandCIJ_bd_cpp(N, Kdir);
bct_test("makerandCIJ_bd", size(CIJ) == [N N] && sum(degrees_dir(CIJ)) == Kdir)

% makerandCIJ_bu
CIJ = makerandCIJ_bu_cpp(N, Kund);
bct_test("makerandCIJ_bu", size(CIJ) == [N N] && sum(degrees_und(CIJ)) == 2 * Kund)

in = floor(rand(N, 1) * floor(N / 3));
out = floor(rand(N, 1) * floor(N / 3));
diff = sum(in) - sum(out);
if diff > 0
	if diff > N
		out += floor(diff / N);
		diff = mod(diff, N);
	end
	out(randperm(N)(1:diff)) += 1;
else
	diff = -diff;
	if diff > N
		in += floor(diff / N);
		diff = mod(diff, N);
	end
	in(randperm(N)(1:diff)) += 1;
end

% makerandCIJdegreesfixed
cij = makerandCIJdegreesfixed_cpp(in, out);
[id od] = degrees_dir(cij);
bct_test("makerandCIJdegreesfixed", in == id' && out == od')

bct_test_teardown
