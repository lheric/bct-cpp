source bct_test_setup.m

subtest = 1;
printf("\nCentrality\n")
source bct_test_centrality.m
printf("\nClustering\n")
source bct_test_clustering.m
printf("\nDegree\n")
source bct_test_degree.m
printf("\nDistance\n")
source bct_test_distance.m
printf("\nModularity\n")
source bct_test_modularity.m
printf("\nMotifs\n")
source bct_test_motifs.m
printf("\nReference\n")
source bct_test_reference.m

printf("\nFailures: %d\n", failures)
clear;
