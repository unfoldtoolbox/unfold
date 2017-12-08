function [spl,nanlist] = dc_designmat_spline(varargin)
% Helper function to generate spline-part of designmatrix
%
% splinecfg{1}  = 'splineName'
% splinecfg{2}  = numberOfSplines
% splinecfg{3}  = 'type', default= bernstein => b-spline
%
%
% cfg.splinespacing
% cfg.codingschema
%
cfg = finputcheck(varargin,...
    {'name',   'string', [], 'spline_default';...
    'paramValues','',[],[];...
    'nsplines','integer',[],[];...
    'knotsequence','real',[],[];...
    'splinespacing','string',{'linear','quantile'},'quantile';
    'type','string',{'default','cyclical','custom'},'default';
    'customfunction','',[],[];
    },'mode','ignore');


assert(~isempty(cfg.nsplines) | ~isempty(cfg.knotsequence),'you need to specify either number of splines or knotsequence')
assert(~isempty(cfg.paramValues),'paramValues were empty')
assert(~all(isnan(cfg.paramValues)),'all paramValues are nans')


spl = [];
spl.paramValues = cfg.paramValues;
spl.nSplines = cfg.nsplines;

%spl.range = range(spl.paramValues);
splmin = min(spl.paramValues);
splmax= max(spl.paramValues);
spl.name = cfg.name;

spl.type = cfg.type;




knots = [];


if isempty(cfg.knotsequence)
    switch cfg.splinespacing
        case 'linear'
            spl.knots =  linspace(splmin,splmax,spl.nSplines-2);
            
        case 'quantile'
            spl.knots =  quantile(spl.paramValues,linspace(0,1,spl.nSplines));
            
        otherwise
            error('wrong cfg.splinespacing. expected linear or quantile')
    end
else
    spl.knots =  knotsequence;
end

if strcmp(spl.type,'default')
    spl.splinefunction = @default_spline;
    
elseif strcmp(spl.type,'cyclical')
    spl.splinefunction = @cyclical_spline;
    
elseif strcmp(spl.type,'custom')
    assert(isa(cfg.customfunction, 'function_handle'),'for custom type one need to define a customfunction')
    spl.splinefunction = cfg.customfunction;
    
else
    error('unknown spline type (should be checked earlier)')
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

nanlist = isnan(spl.paramValues);

spl.X(nanlist,:) = 0; % remove the nans
spl.nSplines = size(spl.X,2);

% give the splines "good" names
% we find the parameter value that maximizes a spline and take this as an
% identifier for the spline. Not perfect, but a good approximation
[~,I] = max(Xspline,[],1);
maxSplVal = spl.paramValues(I);

% round to two significant digits:
tmpSplVal = 2-floor(log10(abs(maxSplVal)));
tmpSplVal(tmpSplVal<0) = 0;
rawColnames = repmat({spl.name},1,spl.nSplines);

spl.colnames = cellfun(@(x,signPoint,y)sprintf('%s_%.*f',x,signPoint,y),rawColnames,num2cell(tmpSplVal)',num2cell(maxSplVal)','UniformOutput', 0);

