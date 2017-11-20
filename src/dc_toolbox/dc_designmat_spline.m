function [spl,nanlist] = dc_designmat_spline(t,splinecfg,cfg)
% Helper function to generate spline-part of designmatrix
%
% splinecfg{1}  = 'splineName'
% splinecfg{2}  = numberOfSplines
%
%
% cfg.splinespacing
% cfg.codingschema
%
spl = [];
spl.paramValues = t{:,splinecfg{1}};
spl.nSplines = splinecfg{2};
%spl.range = range(spl.paramValues);
splmin = min(spl.paramValues);
splmax= max(spl.paramValues);
spl.name = [splinecfg{1}];
% We find the smallest increase to get an idea of the resolution we
% need to define the splines

nanlist = [];

%         spl.spline2val = linspace(spl.min,spl.max,100*spl.nSplines);
%spl.spline2val = spl.min:min(b):spl.max;
for k = 1:length(spl.paramValues)
    if isnan(spl.paramValues(k)) && ~isempty(t.type{k})
        
        
        nanlist = [nanlist k];
        %                 spl.spline2val_idx(k)=1; % this is only temporary, we will remove it after indexing again
        continue
    end
    %             [~,spl.spline2val_idx(k)] = (min(abs(spl.spline2val- spl.paramValues(k))));
end
if ~isempty(nanlist)
    warning('n = %i nan''s found @ %s, event:%i \n currently no support for nan-value in splines. setting it to 0.',length(nanlist),spl.name,k)
end

knots = [];


switch cfg.splinespacing
    case 'linear'
        knots = linspace(splmin,splmax,spl.nSplines-2);
        
    case 'quantile'
        knots = quantile(spl.paramValues,linspace(0,1,spl.nSplines-2));
    otherwise
        error('wrong cfg.splinespacing. expected linear or quantile')
end
% we add 3 knots (because cubic, 4th order splines) in the
% beginning and the end.
knots = [repmat(knots(1),1,3) knots repmat(knots(end),1,3)];
spl.knots = knots;


% This functino always removes either first or last spline. We therefore
% need to recover it by running it twice and concatenating
a = Bernstein(spl.paramValues,spl.knots,[],4,[],0);
b = Bernstein(spl.paramValues,spl.knots,[],4,[],1);

paramValuesSpline = a;
paramValuesSpline(b(:)==1) = 1; 


%%
if strcmp(cfg.codingschema,'effects')
    
    minix = get_min(nanmean(spl.paramValues),spl.paramValues);
    [~,killThisSpline] = max(paramValuesSpline(minix,:));
    
    [~,tmp] = max(paramValuesSpline(:,killThisSpline));
    peakAt = spl.paramValues(tmp);
    fprintf('Due to collinearity, removing the spline for the effect %s has its peak at %f\n',spl.name,peakAt)
    fprintf('This does not mean that the model-intercept represents this value!')
else
    killThisSpline = 1;
end


%         else
%             killThisSpline = 1;
%         end
paramValuesSpline(:,killThisSpline) = [];
spl.removedSplineIdx = killThisSpline;
spl.X= paramValuesSpline;
spl.X(nanlist,:) = 0; % remove the nans
spl.nSplines = size(spl.X,2);

% give the splines "good" names
% we find the parameter value that maximizes a spline and take this as an
% identifier for the spline. Not perfect, but a good approximation
[~,I] = max(paramValuesSpline,[],1);
maxSplVal = spl.paramValues(I);

% round to two significant digits:
tmpSplVal = 2-floor(log10(abs(maxSplVal)));
tmpSplVal(tmpSplVal<0) = 0;
rawColnames = repmat({spl.name},1,spl.nSplines);

spl.colnames = cellfun(@(x,signPoint,y)sprintf('%s_%.*f',x,signPoint,y),rawColnames,num2cell(tmpSplVal)',num2cell(maxSplVal)','UniformOutput', 0);

