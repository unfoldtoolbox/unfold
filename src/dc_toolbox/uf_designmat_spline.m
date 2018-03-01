function [EEG,spl,nanlist] = uf_designmat_spline(EEG,varargin)
% Helper function to generate spline-part of designmatrix
%
%
cfg = finputcheck(varargin,...
    {'name',   'string', [], 'spline_default';...
    'paramValues','',[],[];...
    'nsplines','integer',[],[];...
    'knotsequence','real',[],[];...
    'splinespacing','string',{'linear','quantile'},'quantile';
    'splinefunction','','','default';
    },'mode','ignore');


assert(~isempty(cfg.nsplines) | ~isempty(cfg.knotsequence),'you need to specify either number of splines or knotsequence')
assert(~isempty(cfg.paramValues),'paramValues were empty')
assert(~all(isnan(cfg.paramValues)),'all paramValues are nans')

if size(cfg.paramValues,2) == 1
    cfg.paramValues = cfg.paramValues';
end



spl = [];
spl.paramValues = cfg.paramValues;
spl.nSplines = cfg.nsplines;

%spl.range = range(spl.paramValues);
splmin = min(spl.paramValues);
splmax= max(spl.paramValues);
spl.name = cfg.name;






knots = [];


if isempty(cfg.knotsequence)
    switch cfg.splinespacing
        case 'linear'
            spl.knots =  linspace(splmin,splmax,spl.nSplines-2);
            
        case 'quantile'
            spl.knots =  quantile(spl.paramValues,linspace(0,1,spl.nSplines-2));
            
        otherwise
            error('wrong cfg.splinespacing. expected linear or quantile')
    end
else
    spl.knots =  cfg.knotsequence;
end


if ~ischar(cfg.splinefunction)
    assert(isa(cfg.splinefunction, 'function_handle'),'for custom type one need to define a splinefunction')
    spl.splinefunction = cfg.splinefunction;
elseif strcmp(cfg.splinefunction,'default')
    spl.splinefunction = @default_spline;
    
elseif strcmp(cfg.splinefunction,'cyclical')
    spl.splinefunction = @cyclical_spline;
    
else
    error('unknown spline type')
end

Xspline =  spl.splinefunction(spl.paramValues,spl.knots);

%%
% if strcmp(cfg.codingschema,'effects')

minix = get_min(nanmean(spl.paramValues),spl.paramValues);
[~,killThisSpline] = max(Xspline(minix,:));

[~,tmp] = max(Xspline(:,killThisSpline));
peakAt = spl.paramValues(tmp);
fprintf('The spline that got removed due to collinearity in the basis set (as intended) for the effect %s has its peak at %f\n',spl.name,peakAt)
fprintf('This does not mean that the event-intercept represents this value! \n')
% else
%     killThisSpline = 1;
% end


%         else
%             killThisSpline = 1;
%         end
Xspline(:,killThisSpline) = [];
spl.removedSplineIdx = killThisSpline;
spl.X= Xspline;

% Where the event is not defined a NAN appears, we have to set those to 0


spl.nSplines = size(spl.X,2);

% give the splines "good" names
% we find the parameter value that maximizes a spline and take this as an
% identifier for the spline. Not perfect, but a good approximation
[~,I] = max(Xspline,[],1);
maxSplVal = spl.paramValues(I);

% round to two significant digits:
tmpSplVal = 2-floor(log10(abs(maxSplVal)));
tmpSplVal(isinf(tmpSplVal)) = 0;
tmpSplVal(tmpSplVal<0) = 0;
rawColnames = repmat({spl.name},spl.nSplines,1);

spl.colnames = cellfun(@(x,signPoint,y)sprintf('%s_%.*f',x,signPoint,y),rawColnames,num2cell(tmpSplVal)',num2cell(maxSplVal)','UniformOutput', 0);


%% Add the spline to the EEG-data

EEG.deconv.splines{end+1} = spl;

nanlist = isnan(spl.paramValues);
spl.X(nanlist,:) = 0; % remove nan-entries from splines from designmatrix (for the splines they were removed already)
EEG.deconv.X = [EEG.deconv.X spl.X]; % add spline columns


EEG.deconv.colnames = [EEG.deconv.colnames  spl.colnames'];
EEG.deconv.variablenames = [EEG.deconv.variablenames {spl.name}];


EEG.deconv.cols2variablenames = [EEG.deconv.cols2variablenames repmat(length(EEG.deconv.variablenames),1, spl.nSplines)];
EEG.deconv.cols2eventtypes = [EEG.deconv.cols2eventtypes repmat(EEG.deconv.cols2eventtypes(1),1,size(spl.X,2))];
EEG.deconv.variabletypes = [EEG.deconv.variabletypes 'spline'];
%predType = [predType repmat({'spline'},1,spl.nSplines)];

