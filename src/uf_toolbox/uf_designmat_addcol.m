function EEG = uf_designmat_addcol(EEG,newcol,label,variablettype)
%uf_designmat_addcol
% Adds a single custom column to the unfold-Designmat "Xdc". This is
% sometimes useful to add e.g. continuous predictors manually.
% Note that this is somewhat experimental, and not all functions support
% this.
%
%Arguments:
%   newrow (array): The column(s) to add to the Xdc designmat
%   label(string): The label/identifier of the column
%   eventtype(string): (optional, default nan) the eventtype. E.g. for trf
%   this should be {'trf'}. If you do not manually timeexpand, this should
%   be nan(1). Else the glmfit function does not know to ignore this
%   column while reshaping.
%Return:
%   EEG-Struct
%   * unfold.Xdc added column
%   * unfold.colnames added label


assert(isfield(EEG.unfold,'Xdc'),'could not find unfold.Xdc, run uf_timeexpandDesignmat before')
assert(size(EEG.unfold.Xdc,1) == length(newcol),'New column does not have same size as unfold.Xdc')

if nargin == 3
    variablettype = nan;
else
    assert(ischar(variablettype),'Eventtype has to be a string')

end

% in case we add a whole set of new columns, add them to the same event
if length(size(newcol))==2 && ~any(size(newcol) == 1)
    n = size(newcol,2);
else
    n = 1;
end
% First add a dummy-entry to the fields of "X"
EEG.unfold.X(:,end+1) = nan(size(EEG.unfold.X,1),1);
EEG.unfold.variabletypes(end+1) = {variablettype};
EEG.unfold.variablenames(end+1) = {label};
EEG.unfold.colnames(end+1) = {label};
EEG.unfold.cols2variablenames(end+1) = length(EEG.unfold.variablenames);
EEG.unfold.eventtypes(end+1) = {{nan}};

EEG.unfold.cols2eventtypes(end+1) = length(EEG.unfold.eventtypes);

% Now add the entries to "Xdc"
% Add the new col (cols if timeexpanded set) to Xdc
EEG.unfold.Xdc(:,end+1:(end+n)) = newcol;
% book keeping
EEG.unfold.Xdc_terms2cols(end+1:(end+n)) = length(EEG.unfold.colnames);


if ~iscell(variablettype)
    variablettype = {variablettype};
end

