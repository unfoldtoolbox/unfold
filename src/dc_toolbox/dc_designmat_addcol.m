function EEG = dc_designmat_addcol(EEG,newrow,label)
%dc_designmat_addcol
% Adds a single custom column to the deconv-Designmat "Xdc"
%
%Arguments:
%   newrow (array): The column to add to the Xdc designmat
%   label(string): The label/identifier of the column
%
%Return:
%   EEG-Struct
%   * deconv.Xdc added column
%   * deconv.colnames added label


assert(isfield(EEG.deconv,'Xdc'),'could not find deconv.Xdc, run dc_timeexpandDesignmat before')
assert(size(EEG.deconv.Xdc,1) == length(newrow),'New row does not have same size as deconv.Xdc')

EEG.deconv.Xdc(:,end+1) = newrow;
EEG.deconv.colnames(end+1) = {label};
EEG.deconv.Xdc_terms2cols(end+1) = length(EEG.deconv.colnames);
EEG.deconv.eventtype(end+1) = {nan};
EEG.deconv.cols2eventtype(end+1) = length(EEG.deconv.eventtype);