function [EEG,spl,nanlist] = uf_designmat_spline(EEG,varargin)
%Helper function to generate spline-part of designmatrix
%
%Argument:
%   cfg.name(string): (optional, default: "spline_default") A name for the
%           spline predictor
%
%   cfg.paramValues(double): The values of the predictor on which the
%           splines should be calculated and evaluated.
%           E.g. [-3 4,1,2,3, ... 4]
%
%   cfg.nsplines(integer): number of splines to use. too many lead to
%               overfitting, to few to underfitting. This number will be 
%               transformed into the number of knots later. As different 
%               spline functions have different requiements to the number 
%               of knots, we specify this number and fix the number of
%               knots later.
%
%   cfg.knotsequence(real): optional (if nsplines is specified). Give the
%               sequence of knots explicitly (else they are put on the
%               quantiles or linearly (see splinespacing). An example would
%               be [0,1,2,3,5,10,11,12,13]. This example could make sense
%               if there is lots of data at predictor 0-5 and again at
%               10-13.  To define the knotsequence explicitly is also
%               useful when you directly want to estimate the same betas
%               for all subjects. But beware of subject-specific ranges, not all
%               subjects have the same range in their covariates.
%
%   cfg.splinespacing(string): (quantile) linear or quantile. The spacing
%         of the knots along the
%
%   cfg.splinefunction (functionhandle): You can specify your own spline
%        function. This in principle also allows to make use of polynomial
%        regression
%
%Returns:
% EEG structure
%   * unfold.X: new entries for the spline
%   * unfold.splines: new entrie for the spline
%   * in addition update to unfold: colnames, variablenames,cols2variablenames,cols2eventtypes,variabletypes
% spl: Same as EEG.unfold.spline{end}
% nanlist: the paramValues that were nan (same as 'isnan(spl.paramValues)' )
%
%Example:
%   EEG = dc_designmat_spline(EEG,'name','splineA','paramValues',[EEG.event.splineB],'nsplines',10,'splinespacing','linear');
%   EEG = dc_designmat_spline(EEG,'name','splineB','paramValues',[EEG.event.splineB],'knotsequence',linspace(0,2*pi,15),'splinefunction','cyclical');
%

cfg = finputcheck(varargin,...
    {'name',   'string', [], 'spline_default';...
    'paramValues','',[],[];...
    'nsplines','integer',[],[];...
    'knotsequence','real',[],[];...
    'splinespacing','string',{'linear','quantile'},'quantile';
    'splinefunction','','','default';
    'cyclical_bounds','real','',[];
    },'mode','error');


assert(~isempty(cfg.nsplines) | ~isempty(cfg.knotsequence),'you need to specify either number of splines or knotsequence')
assert(~isempty(cfg.paramValues),'paramValues were empty')

spl = [];
spl.nSplines = cfg.nsplines; 
spl.name = cfg.name;

assert(~all(isnan(cfg.paramValues(:))),'all paramValues are nans')
if strcmp(cfg.splinefunction,'2D')
    if size(cfg.paramValues,2) == 2
        cfg.paramValues = cfg.paramValues';
    end
    splmin = min(cfg.paramValues,[],2);
    splmax= max(cfg.paramValues,[],2);
else
    if size(cfg.paramValues,2) == 1
        cfg.paramValues = cfg.paramValues';
    end
    splmin = min(cfg.paramValues);
    splmax= max(cfg.paramValues);
end
spl.paramValues = cfg.paramValues;


%spl.range = range(spl.paramValues);


% In case of cyclical spline we need to wrap around the spl.paramValues by
% user specified limits before calculating the knot sequence.
if  strcmp(cfg.splinefunction,'cyclical')
    
    if ~isempty(cfg.cyclical_bounds)
        fprintf('Wrapping the parametervalues around the specified lower/upper bound\n')
        
        lbound = cfg.cyclical_bounds(1);
        ubound = cfg.cyclical_bounds(2);
        
        x = spl.paramValues;
        x(x > ubound) = lbound + (x(x > ubound) - ubound); % (ubound - lbound)
        x(x < lbound) = ubound - (lbound - x(x < lbound)); % (ubound - lbound)
        spl.paramValues = x;
        
    elseif isempty(cfg.knotsequence)
        warning('You specified a cyclical predictor, but neither specified a knotsequence nor a lower/upper-bound. We will assume that the lower/upper bound is min/max of your predictor')
        
    end
    if isempty(cfg.knotsequence)
        warning('You specified a cyclical predictor, but did not specify an explicit knotsequence. We will apply the splinespacing function (linear/quantile) on the parameter values to find the knotsequence')
    end
end




if strcmp(cfg.splinefunction,'cyclical')
    %the spline function removes two splines, thus we have to add one of
    %them
    spl.nSplines = spl.nSplines + 1;
end
    

if strcmp(cfg.splinefunction,'default') || strcmp(cfg.splinefunction,'2D')  % jpo 15.04.2018: not sure about this
    % the function adds two splines to the knotsequence, we temporaly
    % remove them 
    spl.nSplines = spl.nSplines-2;
end

knots = [];
if isempty(cfg.knotsequence)
    switch cfg.splinespacing
        case 'linear'
            if strcmp(cfg.splinefunction,'2D')
                spl.knots       = linspace(splmin(1),splmax(1),spl.nSplines);
                spl.knots(2,:)  = linspace(splmin(2),splmax(2),spl.nSplines);
            else
                spl.knots =  linspace(splmin,splmax,spl.nSplines);
            end
        case 'quantile'
            if strcmp(cfg.splinefunction,'2D')
                spl.knots      =  quantile(spl.paramValues(1,:),linspace(0,1,spl.nSplines));
                spl.knots(2,:) =  quantile(spl.paramValues(2,:),linspace(0,1,spl.nSplines));
            else
                spl.knots =  quantile(spl.paramValues,linspace(0,1,spl.nSplines));
            end
        otherwise
            error('wrong cfg.splinespacing. expected linear or quantile')
    end
else
    spl.knots =  cfg.knotsequence;
end

if strcmp(cfg.splinefunction,'default')
    % and we add them again
    spl.nSplines = spl.nSplines+2;
elseif strcmp(cfg.splinefunction,'cyclical')
    spl.nSplines = spl.nSplines-1;
end

if ~ischar(cfg.splinefunction)
    assert(isa(cfg.splinefunction, 'function_handle'),'for custom type one need to define a splinefunction')
    spl.splinefunction = cfg.splinefunction;
    
elseif strcmp(cfg.splinefunction,'default') || strcmp(cfg.splinefunction,'2D')
    spl.splinefunction = @default_spline;
    
elseif strcmp(cfg.splinefunction,'cyclical')
    spl.splinefunction = @cyclical_spline;
    
else
    error('unknown spline type')
end

if strcmp(cfg.splinefunction,'2D')
    Xspline1 =  spl.splinefunction(spl.paramValues(1,:),spl.knots(1,:));
    Xspline2 =  spl.splinefunction(spl.paramValues(2,:),spl.knots(2,:));
    for iD = 1:size(Xspline1,1)
        aux = Xspline1(iD,:)'*Xspline2(iD,:);
        Xspline(iD,:) = aux(:)';
    end
else
    Xspline =  spl.splinefunction(spl.paramValues,spl.knots);
end

%%
% if strcmp(cfg.codingschema,'effects')
if strcmp(cfg.splinefunction,'2D')
    [~,minix] = min(sqrt((nanmean(spl.paramValues(1,:))-spl.paramValues(1,:)).^2+(nanmean(spl.paramValues(2,:))-spl.paramValues(2,:)).^2)); % 
else
    minix = get_min(nanmean(spl.paramValues),spl.paramValues);
end
[~,killThisSpline] = max(Xspline(minix,:));
[~,tmp] = max(Xspline(:,killThisSpline));

if strcmp(cfg.splinefunction,'2D')
    peakAt = spl.paramValues(:,tmp);
    fprintf('The spline that got removed due to collinearity in the basis set (as intended) for the effect %s has its peak at %f,%f\n',spl.name,peakAt(1),peakAt(2))
else
    peakAt = spl.paramValues(tmp);
    fprintf('The spline that got removed due to collinearity in the basis set (as intended) for the effect %s has its peak at %f\n',spl.name,peakAt)
end
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
if strcmp(cfg.splinefunction,'2D')
    maxSplVal = spl.paramValues(:,I);  
else
    maxSplVal = spl.paramValues(I);
end
% round to two significant digits:
tmpSplVal = 2-floor(log10(abs(maxSplVal)));
tmpSplVal(isinf(tmpSplVal)) = 0;
tmpSplVal(tmpSplVal<0) = 0;
rawColnames = repmat({spl.name},spl.nSplines,1);

if strcmp(cfg.splinefunction,'2D')
    spl.colnames = cellfun(@(x,signPoint1,y,signPoint2,z)sprintf('%s_%.*f,%.*f',x,signPoint1,y,signPoint2,z),rawColnames,num2cell(tmpSplVal(1,:))',num2cell(maxSplVal(1,:))',num2cell(tmpSplVal(2,:))',num2cell(maxSplVal(2,:))','UniformOutput', 0);
else
    spl.colnames = cellfun(@(x,signPoint,y)sprintf('%s_%.*f',x,signPoint,y),rawColnames,num2cell(tmpSplVal)',num2cell(maxSplVal)','UniformOutput', 0);
end

%% Add the spline to the EEG-data

EEG.unfold.splines{end+1} = spl;

nanlist = isnan(spl.paramValues);
% spl.X(nanlist,:) = 0; % remove nan-entries from splines from designmatrix (for the splines they were removed already)
EEG.unfold.X = [EEG.unfold.X spl.X]; % add spline columns


EEG.unfold.colnames = [EEG.unfold.colnames  spl.colnames'];
EEG.unfold.variablenames = [EEG.unfold.variablenames {spl.name}];


EEG.unfold.cols2variablenames = [EEG.unfold.cols2variablenames repmat(length(EEG.unfold.variablenames),1, spl.nSplines)];
EEG.unfold.cols2eventtypes = [EEG.unfold.cols2eventtypes repmat(EEG.unfold.cols2eventtypes(1),1,size(spl.X,2))];
EEG.unfold.variabletypes = [EEG.unfold.variabletypes 'spline'];
%predType = [predType repmat({'spline'},1,spl.nSplines)];

