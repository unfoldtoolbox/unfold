function [] = test_addmarginals()

%%

testCase = [15]
    
EEG = simulate_test_case(testCase,'noise',0,'basis','box');
cfgDesign = [];

cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1 + cat(conditionA) + continuousA','y~1'};
cfgDesign.eventtypes = {{'stimulus2'},{'stimulus3'}};

EEG = uf_designmat(EEG,cfgDesign);
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.5 1.5]);

EEG = uf_glmfit(EEG);

ufresult = uf_condense(EEG);

ufresult_average = uf_predictContinuous(ufresult,'auto_method','average');
ufresult_predict = uf_predictContinuous(ufresult,'predictAt',{{'continuousA',[0 1]},{'contB',[0 1]},{'splineA',[-1 0 1]}});
ufresult_marginal = uf_addmarginal(ufresult_predict);

% the stimulus3 sould be unaffected, it's just an intercept
assert(1==near(ufresult_marginal.beta(:,:,5),ufresult_predict.beta(:,:,5))) % this column should not have been impacted.

% the stimulus 2, parameter 3 is continuous = 0, so it should be the same
% as the intercept
assert(1==near(ufresult_marginal.beta(:,:,3), ufresult.beta(:,:,1)))

% the stimulus 2, parameter 1 should be at continuous = mean(continuous)
% {that's 50} so it should be continuous-average + intercept
assert(1==near(ufresult_marginal.beta(:,:,1), sum(ufresult_average.beta(:,:,[1,3]),3)))



end

