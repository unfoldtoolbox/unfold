function EEG = dc_designmat_addcol(EEG,newrow,label)
%dc_designmat_addcol
% Adds a single custom column to the deconv-Designmat "dcX"
%
%Arguments:
%   newrow (array): The column to add to the dcX designmat
%   label(string): The label/identifier of the column
%
%Return:
%   EEG-Struct
%   * deconv.dcX added column
%   * deconv.colnames added label


assert(isfield(EEG.deconv,'dcX'),'could not find deconv.dcX, run dc_timeexpandDesignmat before')
assert(size(EEG.deconv.dcX,1) == length(newrow),'New row does not have same size as deconv.dcX')

EEG.deconv.dcX(:,end+1) = newrow;
EEG.deconv.colnames(end+1) = {label};
EEG.deconv.dcX_termidx(end+1) = length(EEG.deconv.colnames);
EEG.deconv.eventtype(end+1) = {nan};
EEG.deconv.col2eventtype(end+1) = length(EEG.deconv.eventtype);