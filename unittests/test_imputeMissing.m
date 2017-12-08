function test_imputeMissing()
EEG = simulate_test_case(5,'noise',0,'basis','box');
for e = 10:20
    EEG.event(e).continuousA = nan(1);
end
EEG = dc_designmat(EEG,'formula','y~1+continuousA','eventtypes','stimulusA');


check_nan(dc_imputeMissing(EEG,'method','mean'))
check_nan(dc_imputeMissing(EEG,'method','median'))
check_nan(dc_imputeMissing(EEG,'method','marginal'))
check_nan(dc_imputeMissing(EEG,'method','drop'))
EEG2 = dc_imputeMissing(EEG,'method','drop');


%% Test that glmfit crashes when nans are there

try
    EEG = dc_timeexpandDesignmat(EEG,'timelimits',[0 1]);
    error('this should have thrown an error')
catch
    EEG2 = dc_timeexpandDesignmat(EEG2,'timelimits',[0 1]);
end

try
    dc_glmfit(EEG2,'method','lsmr');
    error('this should have thrown an error')
catch
end

display('imputation methods tests without error')
end


function check_nan(EEG)
if any(isnan(EEG.deconv.X(:)))
    error('imputeMissing has some kind of problem')
end
end