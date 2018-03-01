% How to interprete results
% - all odd numbers do not have overlap, thus the normal ERP and the
%   deconvolved one should be identical and the same to the original signal
%
%
%Testing one event
% 1 	- One Intercept only no overlap
% 2 	- One Intercept only with overlap
% 3 	- One 1x2 Factor with an effect but not on timing
% 4 	- One 1x2 Factor with no effect except on timing
% 5 	- One continuous Factor with abs(slope) > 0 but not on timing
% 6 	- One continuous Factor with 0-slope but slope on timing
%
% b) Testing splines
% 7 	- intercept plus simple cubic effect, no overlap
% 8 	- intercept plus simple cubic effect, with overlap
%
%
% c) Testing multiple events
% 9    	- Two events, no overlap
% 10   	- Two events, with overlap
% 11  	- Three events, no overlap
% 12   	- Three events, with overlap
% 13   	- Three events, one intercept, one 1x2 and one continuos, no overlap
% 14   	- Three events, one intercept, one 1x2 and one continuos, with overlap

testCase = 15


EEG = simulate_test_case(testCase,'noise',0,'basis','hanning');

cfgDesign = [];
cfgDesign.eventtypes = {'stimulusA'};
cfgDesign.codingschema = 'reference';
switch testCase
    case {1,2}
        cfgDesign.formula = 'y ~ 1';
    case {3,4}
        cfgDesign.formula = 'y ~ 1 + conditionA';
        cfgDesign.categorical = {'conditionA'};
    case {5,6}
        cfgDesign.formula = 'y ~ 1+ continuousA';
    case {7,8}
        cfgDesign.formula = 'y ~ 1';
        cfgDesign.spline = {{'splineA',10}}; % In addition use one spline
    case {9,10}
        cfgDesign.formula = {'y~1','y~1'};
        cfgDesign.eventtypes = {'stimulusA','stimulusB'};
    case {11,12}
        cfgDesign.formula = {'y~1','y~1','y~1'};
        cfgDesign.eventtypes = {'stimulusA','stimulusB','stimulusC'};
    case {13,14}
        cfgDesign.formula = {'y~1','y~1+conditionA','y~1+continuousA'};
        cfgDesign.eventtypes = {'stimulusA','stimulusB','stimulusC'};
        cfgDesign.categorical = {'conditionA'};
    case {15}
        cfgDesign.formula   = {'y~1',       'y~1+cat(conditionA)*continuousA', 'y~1+spl(splineA,5)+spl(splineB,5)+continuousA'};
        cfgDesign.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};
end

EEG = uf_designmat(EEG,cfgDesign);
% error
timelimits = [-1.5,2.5];
EEG = uf_timeexpandDesignmat(EEG,'timelimits',timelimits,'method','full','timeexpandparam',7);
EEG= uf_glmfit(EEG,'method','lsmr');
% EEG2= uf_glmfit(EEG,'method','glmnet');

%%
EEG = uf_epoch(EEG,'timelimits',timelimits);
EEG = uf_glmfit_nodc(EEG); %does not overwrite

%%
unfold = uf_condense(EEG);

multWith = ones(1,size(EEG.deconv.X,2));
for col = 1:size(EEG.deconv.X,2)
    ix = ismember({EEG.urevent.type},EEG.deconv.eventtypes{EEG.deconv.cols2eventtypes(col)});
    multWith(col) = mean(EEG.deconv.X(ix,col),1);
end
%%

figure
subplot(2,1,1)
plot(unfold.times,bsxfun(@times,squeeze(unfold.beta),multWith),'-x'),hold all
% plot(unfold.times,squeeze(unfold.beta),'-x'),hold all
plot(EEG.sim.sig.time,EEG.sim.separateSignal','-ok')
title('unfold vs. orig')

subplot(2,1,2)
plot(unfold.times,bsxfun(@times,squeeze(unfold.beta_nodc),multWith),'-x'),hold all
% plot(unfold.times,squeeze(unfold.beta_nodc),'-x'),hold all
plot(EEG.sim.sig.time,EEG.sim.separateSignal','-ok')
title('epoched vs. orig')
%% draw splinethings
unfold = uf_condense(EEG);

cfg = [];
cfg.auto_method = 'linear'; %default quantile
cfg.auto_n = 20;
cfg.convertSplines = 1;

unfold = uf_getParam(unfold,cfg);


cfg = [];
cfg.channel = 1;
cfg.sameyaxis = 'all';
cfg.deconv = -1;
cfg.plotSeparate = 'event';
cfg.plotParam = {'3_(Intercept)','3_continuousA','splineA','splineB'};
cfg.add_marginal = 0;
ax = uf_plotParam(unfold,cfg);
