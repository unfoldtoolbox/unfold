function [splIdxListAll,paramList] = uf_getSplineidx(EEG)
%returns a list of the splines in EEG and a parameter list where only the
% first index of the spline remains
%
%Arguments:
%   EEG: eeglab set with EEG.deconv.spl.
%
%Returns:
% splIdxList:  Empty if no splines, else contains the column indices of the designmatrix X where there are splines
% paramList: Returns a list with the number of parameters (The number of columns in X), but only contains a single one for the spline.
%
%*Example:*
% A designmatrix with 2 normal predictors, a spline with 5 knots
%| and another normal predictor.
%| Call [splIdxList,paramList] = uf_get_splines(EEG)
%| splIdxList = [3 4 5 6 7]
%| paramList= [1 2 3 8]


assert(isfield(EEG,'deconv')&&isfield(EEG.deconv,'splines'),'could not find EEG.deconv.splines')

numParam = size(EEG.deconv.X,2);
splIdxListAll = [];
paramList = 1:numParam;


splineVar = find(strcmp(EEG.deconv.variabletypes,'spline'));
splineCol = EEG.deconv.cols2variablenames;
splineCol(~ismember(EEG.deconv.cols2variablenames,splineVar)) = 0;

if isfield(EEG.deconv,'splines') && ~isempty(EEG.deconv.splines)
    for splIdx = 1:length(splineVar)
        splIdxList = find(splineCol == splineVar(splIdx));
%         splIdxList = find(cellfun(@check_spline,EEG.deconv.colnames(:),repmat({splIdx},length(EEG.deconv.colnames),1)));
        paramList = setdiff(paramList,splIdxList(2:end)); % we want to evaluate the sopline signal only once
        splIdxListAll = [splIdxListAll splIdxList];
    end
else
    splIdxListAll = [];
end


end
