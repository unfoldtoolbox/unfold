function test_splines()
% Generate data with only a single timepoint
% apply a non-linear function and try to recover it.
%%
cfgSim = [];
cfgSim.plot = 1;
for type = {'default','cyclical','custom','cyclical_formula','2D'}
    
    cfgSim.type = type{1};
    
    
    
    %% Generate data
    switch cfgSim.type
        case {'cyclical','cyclical_formula'}
            spl.function = @cyclical_spline;
            spl.values = linspace(0,4*pi,100);
            datafunction = @(x)(sin(x+pi/2)+0.5*sin(3*x)); % something where we need phase diff :)
        case {'default','custom'}
            
            spl.function = @default_spline;
            spl.values = linspace(-10,10,100);
            datafunction = @(x)x.^2; % something where we need phase diff :)
        case {'2D'}
            spl.function = []; % autoselect
            [x,y]= meshgrid(linspace(-10,10,100),linspace(-10,10,100));
            spl.values  = [x(:) y(:)];
            datafunction = @(x)(x(:,1).^2 + 0.05*x(:,2).^3)';
    end
    % Xspline = splinefunction(splinevalues,[linspace(0,2*pi,10)]);
    
  EEG = test_spline_genEEG(spl,datafunction,cfgSim.type);
    
    %% fit
    switch cfgSim.type
        case 'default'
            EEG = uf_designmat(EEG,'eventtypes','stimulus','formula','y~1+spl(splineA,10)');
        case 'cyclical'
            EEG = uf_designmat(EEG,'eventtypes','stimulus','formula','y~1');
            EEG = uf_designmat_spline(EEG,'name','splineA','paramValues',[EEG.event.splineA],'knotsequence',linspace(0,2*pi,15),'splinefunction','cyclical');
        case 'custom'
            EEG = uf_designmat(EEG,'eventtypes','stimulus','formula','y~1');
            EEG = uf_designmat_spline(EEG,'name','splineA','paramValues',[EEG.event.splineA],'knotsequence',linspace(-10,10,5),'splinefunction',spl.function);
        case 'cyclical_formula'
            EEG = uf_designmat(EEG,'eventtypes','stimulus','formula','y~1+circspl(splineA,15,0,2*pi)');
        
        case '2D'
            EEG = uf_designmat(EEG,'eventtypes','stimulus','formula','y~1+2dspl(splineA,splineB,4)');
    end
    
    
    EEGepoch = uf_epoch(EEG,'timelimits',[0 1]);
    EEGepoch = uf_glmfit_nodc(EEGepoch);
    ufresult = uf_condense(EEGepoch);
    if strcmp(cfgSim.type,'2D')
        ufresultconverted = uf_predictContinuous(ufresult,'predictAt',{{'splineAsplineB',linspace(-10,10,100),linspace(-10,10,100)}});
    else
        ufresultconverted = uf_predictContinuous(ufresult,'predictAt',{{'splineA',spl.values}});
    end
    dc  =  ufresultconverted.beta_nodc(:,:,1);
    result = squeeze(ufresultconverted.beta_nodc(:,:,2:end) +dc);
    
    if cfgSim.plot
        
        figure('name',cfgSim.type)
        if strcmp(cfgSim.type,'2D')
            %%
%             figure,
            modelled = reshape(result,100,100)';
            data = reshape(EEG.data,100,100);
            subplot(3,1,1)
            imagesc(modelled)
            title('model')
            colorbar
            caxis([-50 150])
            
            subplot(3,1,2)
            imagesc(data)
            caxis([-50 150])
            colorbar
            title('data')
            
            subplot(3,1,3)
            imagesc(modelled-data)
            title('residuals = model - data')
            caxis([-20 20])
        
            colorbar
        else
        subplot(2,1,1)
        Xspline = spl.function(spl.values,EEG.unfold.splines{1}.knots);
        Xspline(:,EEG.unfold.splines{1}.removedSplineIdx) = [];
        plot(dc + Xspline*squeeze(ufresult.beta_nodc(:,:,2:end)),'Linewidth',2)
        hold on
        plot(dc + Xspline.*squeeze(ufresult.beta_nodc(:,:,2:end))',':','Linewidth',2)
        
        % original data vs recovery
        subplot(2,1,2)
        plot(dc+ Xspline*squeeze(ufresult.beta_nodc(:,:,2:end)),'Linewidth',2)
        hold on
        plot(EEG.data,'--','Linewidth',2)
        end
        
    end
    if strcmp(cfgSim.type,'2D')
        assert(sum((modelled(:)-data(:)).^2) < 0.002,'could not recover function!')
    else
        
    assert(sum((EEG.data- result').^2) < 0.001,'could not recover function!')
    end
    
end    
%% Spline & Imputation
% Ticket #46
EEGtmp = simulate_test_case(7,'noise',0,'basis','box');
EEGtmp.event(1).splineA = nan;
EEGtmp = uf_designmat(EEGtmp,'eventtypes','stimulusA','formula','y~1+spl(splineA,4)');
assert(all(isnan(EEGtmp.unfold.X(1,2:end))))

%% 2D Spline & addmarginal
% Ticket #116
spl.function = []; % autoselect
[x,y]= meshgrid(linspace(-10,10,100),linspace(-10,10,100));
spl.values  = [x(:) y(:)];
datafunction = @(x)(x(:,1).^2 + 0.05*x(:,2).^3)';

    
EEGtmp = test_spline_genEEG(spl,datafunction,'2D');
  
EEGtmp = uf_designmat(EEGtmp,'eventtypes','stimulus','formula','y~1+2dspl(splineA,splineB,4)');

EEGtmp = uf_epoch(EEGtmp,'timelimits',[0 1]);
EEGtmp = uf_glmfit_nodc(EEGtmp);

ufresult = uf_condense(EEGtmp);
ufresultconverted = uf_predictContinuous(ufresult,'predictAt',{{'splineAsplineB',linspace(-10,10,100),linspace(-10,10,100)}});
uf_predictContinuous(ufresult)
uf_addmarginal(ufresultconverted)

end
