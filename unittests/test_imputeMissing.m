function test_imputeMissing()
EEG = simulate_test_case(13,'noise',0,'basis','box');
stimCix = find(strcmpi({EEG.event.type},'stimulusC'));
% 3 stimuliC need to be imputed

for e = 5:8
    EEG.event(stimCix(e)).continuousA = nan(1);
end
EEG = dc_designmat(EEG,'formula',{'y~1','y~1+continuousA'},'eventtypes',{{'stimulusA'},{'stimulusC'}});


for type = {'mean','median','marginal','drop'}
    EEGimp = dc_imputeMissing(EEG,'method',type{1});

    switch type{1}
        case 'mean'
            assert(all(nanmean(EEG.deconv.X) - nanmean(EEGimp.deconv.X)< 10^-10),'error, mean imputation showed larger error')
        case 'drop'
            tmp = EEGimp.deconv.X(stimCix(5:8),:);
            assert(all(tmp(:)==0),'error, drop imputation did something wrong')
        otherwise
            
    end
    
check_nan(EEGimp)
end

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