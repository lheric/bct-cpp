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

bct_test_teardown
