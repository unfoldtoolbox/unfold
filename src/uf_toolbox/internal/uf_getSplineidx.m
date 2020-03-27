function [splIdxListAll,paramList] = uf_getSplineidx(EEG)
%returns a list of the splines in EEG and a parameter list where only the
% first index of the spline remains
%
%Arguments:
%   EEG: eeglab set with EEG.unfold.spl.
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


assert(isfield(EEG,'unfold')&&isfield(EEG.unfold,'splines'),'could not find EEG.unfold.splines')

% numParam = size(EEG.unfold.X,2);
paramList= find(cellfun(@(x)~all(isnan(x)),EEG.unfold.variabletypes));
paramList = find(ismember(EEG.unfold.cols2variablenames,paramList));

splIdxListAll = [];


splineVar = find(strcmp(EEG.unfold.variabletypes,'spline'));
splineCol = EEG.unfold.cols2variablenames;
splineCol(~ismember(EEG.unfold.cols2variablenames,splineVar)) = 0;

if isfield(EEG.unfold,'splines') && ~isempty(EEG.unfold.splines)
    for splIdx = 1:length(splineVar)
        splIdxList = find(splineCol == splineVar(splIdx));
%         splIdxList = find(cellfun(@check_spline,EEG.unfold.colnames(:),repmat({splIdx},length(EEG.unfold.colnames),1)));
        paramList = setdiff(paramList,splIdxList(2:end)); % we want to evaluate the sopline signal only once
        splIdxListAll = [splIdxListAll splIdxList];
    end
else
    splIdxListAll = [];
end


end
