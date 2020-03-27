function test_imputeMissing()
EEG = simulate_test_case(13,'noise',0,'basis','box');
stimCix = find(strcmpi({EEG.event.type},'stimulusC'));
% 3 stimuliC need to be imputed

for e = 5:8
    EEG.event(stimCix(e)).continuousA = nan(1);
end
EEG = uf_designmat(EEG,'formula',{'y~1','y~1+continuousA'},'eventtypes',{{'stimulusA'},{'stimulusC'}});


for type = {'mean','median','marginal','drop'}
    EEGimp = uf_imputeMissing(EEG,'method',type{1});

    switch type{1}
        case 'mean'
            assert(all(nanmean(EEG.unfold.X) - nanmean(EEGimp.unfold.X)< 10^-10),'error, mean imputation showed larger error')
        case 'drop'
            tmp = EEGimp.unfold.X(stimCix(5:8),:);
            assert(all(tmp(:)==0),'error, drop imputation did something wrong')
        otherwise
            
    end
    
check_nan(EEGimp)
end

EEG2 = uf_imputeMissing(EEG,'method','drop');


%% Test that glmfit crashes when nans are there

try
    EEG = uf_timeexpandDesignmat(EEG,'timelimits',[0 1]);
    error('this should have thrown an error')
catch
    EEG2 = uf_timeexpandDesignmat(EEG2,'timelimits',[0 1]);
end

try
    uf_glmfit(EEG2,'method','lsmr');
    error('this should have thrown an error')
catch
end

display('imputation methods tests without error')

%% check spline + continuous nan
EEG = simulate_test_case(13,'noise',0,'basis','box');
stimCix = find(strcmpi({EEG.event.type},'stimulusC'));
% 3 stimuliC need to be imputed

for e = 5:8
    EEG.event(stimCix(e)).continuousA = nan(1);
end
EEG = uf_designmat(EEG,'formula',{'y~1','y~1+spl(continuousA,5)'},'eventtypes',{{'stimulusA'},{'stimulusC'}});

assert(all(find(any(isnan(EEG.unfold.X),2)) == stimCix(5:8)'),'splines do not have nans')

%% bug #82
% problem if multiple events are combined like e.g. this: {{'A','B'},'C'}
% with {'y~1','y~1}
EEG = simulate_test_case(13,'noise',0,'basis','box');

stimAix = find(strcmpi({EEG.event.type},'stimulusA'));

for e = 5:8
    EEG.event(stimAix(e)).continuousA = nan(1);
end
EEG = uf_designmat(EEG,'formula',{'y~1+continuousA','y~1+continuousA'},'eventtypes',{{'stimulusA','stimulusB'},{'stimulusC'}});
EEG = uf_imputeMissing(EEG,'method','drop');

end


function check_nan(EEG)
if any(isnan(EEG.unfold.X(:)))
    error('imputeMissing has some kind of problem')
end
end