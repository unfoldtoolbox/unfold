function test_imputeMissing()
EEG = simulate_test_case(5,'noise',0,'basis','box');
for e = 10:20
    EEG.event(e).continuousA = nan(1);
end
EEG = dc_designmat(EEG,'formula','y~1+continuousA','eventtype','stimulusA');


check_nan(dc_imputeMissing(EEG,'method','mean'))
check_nan(dc_imputeMissing(EEG,'method','median'))
check_nan(dc_imputeMissing(EEG,'method','marginal'))
check_nan(dc_imputeMissing(EEG,'method','drop'))

display('imputation methods tests without error')

end


function check_nan(EEG)
if any(isnan(EEG.deconv.X(:)))
    error('imputeMissing has some kind of problem')
end
end