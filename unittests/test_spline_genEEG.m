function EEG = test_spline_genEEG(spl,datafunction,type)
  EEG = eeg_emptyset();
    EEG.srate =1;
    EEG.xmin = 1;
    EEG.xmax = size(spl.values,2);
    EEG.times = EEG.xmin:EEG.xmax * 1000;
    
    EEG.data = zeros(1,length(EEG.times));
    
    EEG.data = datafunction(spl.values);
%     EEG.data = EEG.data + randn(size(EEG.data)).*1;
    
    EEG = eeg_checkset(EEG);
    
    for e = 1:size(EEG.data,2)
        evt = struct();
        evt.latency = e;
        evt.type = 'stimulus';
        
        if strcmp(type,'2D')
            evt.splineA = spl.values(e,1);
            evt.splineB = spl.values(e,2);
        else
            evt.splineA = spl.values(e);
        end
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
end