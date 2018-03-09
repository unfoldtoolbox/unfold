function EEG = uf_designmat_addcol(EEG,newrow,label)
%uf_designmat_addcol
% Adds a single custom column to the unfold-Designmat "Xdc". This is
% sometimes useful to add e.g. continuous predictors manually.
%
%Arguments:
%   newrow (array): The column to add to the Xdc designmat
%   label(string): The label/identifier of the column
%
%Return:
%   EEG-Struct
%   * unfold.Xdc added column
%   * unfold.colnames added label


assert(isfield(EEG.unfold,'Xdc'),'could not find unfold.Xdc, run uf_timeexpandDesignmat before')
assert(size(EEG.unfold.Xdc,1) == length(newrow),'New row does not have same size as unfold.Xdc')

EEG.unfold.Xdc(:,end+1) = newrow;
EEG.unfold.colnames(end+1) = {label};
EEG.unfold.Xdc_terms2cols(end+1) = length(EEG.unfold.colnames);
EEG.unfold.eventtypes(end+1) = {nan};
EEG.unfold.cols2eventtypes(end+1) = length(EEG.unfold.eventtypes);