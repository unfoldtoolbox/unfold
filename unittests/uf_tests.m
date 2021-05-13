function uf_tests()
%This function runs the main dc-functions with a lot of parameter pairings
%(~10.000) and tests whether errors occur
% It also runs the tests in this folder 
% It also checks how large the difference between estimated signal and
% original signal is and throws an error if the difference is too large

%separate tests
test_addmarginals
test_continuousArtifact
test_designmat
test_glmfit
test_imputeMissing
test_splines
test_timeexpandDesignmat
test_timeexpandDesignmat_addTRF
test_checkmodelfit
test_erpimage

end
