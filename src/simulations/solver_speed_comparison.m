

%% The function reseats the random generator (and sets it back)
EEGsim = simulate_test_case(15,'noise',0,'basis','box','datalength',10*60,'srate',500);
EEGsim.data = EEGsim.data  + randn(size(EEGsim.data));
EEGsim.data = repmat(EEGsim.data,50,1);
cfg = [];
cfg.formula   = {'y~1',       'y~1+cat(conditionA)*continuousA', 'y~1+spl(splineA,5)+spl(splineB,5)+continuousA'};
cfg.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};

EEGsim = dc_designmat(EEGsim,cfg);
EEGsim = dc_timeexpandDesignmat(EEGsim,'timelimits',[-0.5 0.5],'method','full');

%%

tAll = table();
for method = {'lsmr','lsqr','matlab'}
   t1 = tic;
   
   EEGsolved = dc_glmfit(EEGsim,'method',method{1},'channel',[1 2]);
   t2 = toc(t1);
   
   t = table(t2,method(1));
   tAll = [tAll;t];
end


%% CPU VS GPU, LSQR vs LSMR

