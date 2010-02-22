bct_test_setup

N = floor(91 * rand()) + 10;
Kdir = floor((N * N - N + 1) * rand());
Kund = floor(Kdir / 2);

% makerandCIJ_dir
CIJ = makerandCIJ_dir(N, Kdir);
bct_test("makerandCIJ_dir", size(CIJ) == [N N] && sum(degrees_dir(CIJ)) == Kdir)

% makerandCIJ_und
CIJ = makerandCIJ_und(N, Kund);
bct_test("makerandCIJ_und", size(CIJ) == [N N] && sum(degrees_und(CIJ)) == 2 * Kund)

bct_test_teardown
