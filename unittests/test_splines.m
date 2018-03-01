function test_splines()
% Generate data with only a single timepoint
% apply a non-linear function and try to recover it.
%%
cfgSim = [];
cfgSim.plot = 1;
for type = {'default','cyclical','custom'}
    
    cfgSim.type = type{1};
    
    
    
    %% Generate data
    switch cfgSim.type
        case 'cyclical'
            spl.function = @cyclical_spline;
            spl.values = linspace(0,4*pi,100);
            datafunction = @(x)(sin(x+pi/2)+0.5*sin(3*x)); % something where we need phase diff :)
        case {'default','custom'}
            
            spl.function = @default_spline;
            spl.values = linspace(-10,10,100);
            datafunction = @(x)x.^2; % something where we need phase diff :)
            
    end
    % Xspline = splinefunction(splinevalues,[linspace(0,2*pi,10)]);
    
    EEG = eeg_emptyset();
    EEG.srate =1;
    EEG.xmin = 1;
    EEG.xmax = size(spl.values,2);
    EEG.times = EEG.xmin:EEG.xmax * 1000;
    
    EEG.data = zeros(1,length(EEG.times));
    
    EEG.data = datafunction(spl.values);
    
    EEG = eeg_checkset(EEG);
    
    for e = 1:size(EEG.data,2)
        evt = struct();
        evt.latency = e;
        evt.type = 'stimulus';
        evt.splineA = spl.values(e);
        
        % test the robustness to NAN values in splines
        if e == 60 && e ==70
            evt.splineA = nan;
        end
        if isempty(EEG.event)
            EEG.event = evt;
        else
            EEG.event(e) = evt ;
        end
    end
    
    %% fit
    switch cfgSim.type
        case 'default'
            EEG = dc_designmat(EEG,'eventtypes','stimulus','formula','y~1+spl(splineA,10)');
        case 'cyclical'
            EEG = dc_designmat(EEG,'eventtypes','stimulus','formula','y~1');
            EEG = dc_designmat_spline(EEG,'name','splineA','paramValues',[EEG.event.splineA],'knotsequence',linspace(0,2*pi,15),'splinefunction','cyclical');
        case 'custom'
            EEG = dc_designmat(EEG,'eventtypes','stimulus','formula','y~1');
            EEG = dc_designmat_spline(EEG,'name','splineA','paramValues',[EEG.event.splineA],'knotsequence',linspace(-10,10,5),'splinefunction',spl.function);
    end
    
    EEGepoch = dc_epoch(EEG,'timelimits',[0 1]);
    EEGepoch = dc_glmfit_nodc(EEGepoch);
    unfold = dc_condense(EEGepoch);
    unfoldconverted = dc_getParam(unfold,'predictAt',{{'splineA',spl.values}});
    dc  =  unfoldconverted.beta_nodc(:,:,1);
    result = squeeze(unfoldconverted.beta_nodc(:,:,2:end) +dc);
    
    if cfgSim.plot
        
        figure('name',cfgSim.type)
        subplot(2,1,1)
        Xspline = spl.function(spl.values,EEG.deconv.splines{1}.knots);
        Xspline(:,EEG.deconv.splines{1}.removedSplineIdx) = [];
        plot(dc + Xspline*squeeze(unfold.beta_nodc(:,:,2:end)),'Linewidth',2)
        hold on
        plot(dc + Xspline.*squeeze(unfold.beta_nodc(:,:,2:end))',':','Linewidth',2)
        
        % original data vs recovery
        subplot(2,1,2)
        plot(dc+ Xspline*squeeze(unfold.beta_nodc(:,:,2:end)),'Linewidth',2)
        hold on
        plot(EEG.data,'--','Linewidth',2)
        
        
    end
    assert(sum((EEG.data - result').^2) < 0.001,'could not recover function!')
end