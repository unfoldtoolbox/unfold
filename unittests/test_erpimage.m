function test_erpimage
testCase = [15];

EEGsim = simulate_test_case(testCase,'noise',1,'basis','posneg','srate',50,'datalength',10*60);
%EEGsim.event({EEGsim.event.type} == "stimulusB") = []
EEGsim.data(2,:) = EEGsim.data(1,:);
cfgDesign = [];


cfgDesign.formula   = {'y~1',       'y~1'};
cfgDesign.eventtypes = {'stimulusA', 'stimulusB'};

cfgDesign.codingschema = 'effects';
cfgDesign.formula   = {'y~1',       'y~1+cat(conditionA)+continuousA', 'y~1'};
cfgDesign.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};

uf = uf_designmat(EEGsim,cfgDesign);
uf = uf_timeexpandDesignmat(uf,'timelimits',[-.3 5.5]);
uf = uf_glmfit(uf);



try
    uf_erpimage(uf) % check for required fields
    error('uf_erpimage should have thrown an error')
catch
end

try
    uf_erpimage(uf,'channel',1,'type','raw','addResiduals',1) % check for required fields
    error('uf_erpimage should have thrown an error')
catch
end
%% default call

uf_erpimage(uf,'channel',[1,2]);
%% Check types
figure
subplot(3,1,1),uf_erpimage(uf,'channel',[1,2],'type','residual','figure',0),title('residual')
subplot(3,1,2),uf_erpimage(uf,'channel',[1,2],'type','modelled','figure',0),title('modelled')
subplot(3,1,3),uf_erpimage(uf,'channel',[1,2],'type','raw','figure',0),title('raw')
%%
figure
subplot(3,1,1),uf_erpimage(uf,'channel',1,'overlap',0,'addResiduals',0,'figure',0),title('addResiduals = 0')
subplot(3,1,2),uf_erpimage(uf,'channel',1,'overlap',0,'addResiduals',1,'figure',0),title('addResiduals = 1')
subplot(3,1,3),uf_erpimage(uf,'channel',1,'overlap',0,'addResiduals',2,'figure',0),title('addResiduals = 2')

%%
figure
plotList = [1,6,7,8,2,10,3, 11,4,12];
cfgPlot = [];
cfgPlot.channel=[1,2];
cfgPlot.figure = 0;
cfgPlot.alignto = 'stimulus2';
cfgPlot.sort_alignto = {'stimulus2'};
cfgPlot.sort_direction = 'forward';
cfgPlot.caxis = [-6,6];
cfgPlot.sort_time = [0 0];
cfgPlot.sort_by = 'continuousA';
cfgPlot.keep = {'continuousA'};


cfgPlot.overlap = 1;
uf_erpimage(uf,cfgPlot)


%% check for residual correctness
figure
plotList = [1,6,7,8,2,10,3, 11,4,12];
cfgPlot = [];
cfgPlot.figure = 1;
cfgPlot.channel=1;
cfgPlot.alignto = 'stimulus2';
cfgPlot.sort_alignto = {'stimulus3'};
cfgPlot.sort_direction = 'forward';
cfgPlot.remove = {'3_(Intercept)'};

cfgPlot.caxis = [-6,6];



subplot(3,4,plotList(1))
cfgPlot.type = 'raw';
uf_erpimage(uf,cfgPlot)
title('raw')

subplot(3,4,plotList(2))
cfgPlot.overlap = 0;
cfgPlot.type = 'residual';
cfgPlot.datafield = 'beta_nodc';

uf_erpimage(uf,cfgPlot)
title('residual overlap = 0, beta_nodc')

subplot(3,4,plotList(3))
cfgPlot.overlap = 0;
cfgPlot.type = 'residual';
cfgPlot.datafield = 'beta';
uf_erpimage(uf,cfgPlot)
title('residual overlap = 0, beta')

subplot(3,4,plotList(4))
cfgPlot.overlap = 1;
cfgPlot.type = 'residual';
cfgPlot.datafield = 'beta';
uf_erpimage(uf,cfgPlot)
title('residual overlap = 1, beta')
k = 5;

for field = {'beta_nodc','beta'}
    for overlap = [0 1]
        for residuals = [0 1]
            if overlap == 1 && field == "beta_nodc"
                continue
            end
            subplot(3,4,plotList(k))
            k = k + 1;
            cfgPlot.type ='modelled';
            
            cfgPlot.overlap= overlap;
            cfgPlot.addResiduals = residuals;
            cfgPlot.datafield = field{1};
            
            uf_erpimage(uf,cfgPlot)
            
            title(sprintf('overlap:%i, residuals:%i, field:%s',overlap,residuals,field{1}))
        end
    end
end