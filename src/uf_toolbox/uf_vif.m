function vif = uf_vif(X)
% The VIF is the diagonal of the inverse of the correlation between
% coefficients.
% This Algorithm is taken from Belsley, D. A., E. Kuh, and R. E. Welsch. Regression Diagnostics. Hoboken, NJ: John Wiley & Sons, 1980.
%
% Input: A designmatrix X with Columns => Predictors, Rows =>
% Trials/Timepoints. Both Xdc and X can be used.
%
% Output: Variance Inflation Factor. Keep in mind that if you use splines,
% a high VIF is to be expected and usually not a problem (as long as the
% solver converges)
%
% Memory: This function needs 3 x N² bytes Memory. So for a 40.000 predictor
% designmatrix, it takes ~4GB of memory
%
% Example: 
% > vif=uf_vif(EEG.unfold.Xdc)
% > plot(vif,'-o')
% Algorithm to calculate initial 'corrx' taken from Matthew Gunn
% https://stackoverflow.com/questions/32106370/corr-with-sparse-matrix-matlab/33904776#33904776

if size(X,2) > 30000
    warning('This function needs 3*n^2 bytes ram, so ~ %.1fGB for your matrix \n',size(X,2).^2*3/1024^3)
end
[n, ~]  = size(X);
Exxprim = full(X'*X)/n; %I'm shocked if this isn't full so let's drop sparse now 
Ex   = full(mean(X))'; %same deal
covx = (Exxprim - Ex*Ex');
stdx = sqrt(diag(covx));
corrx = covx ./ (stdx * stdx');
clear COVX 
clear Exxprim

fprintf('Calculating Inverse of Correlation Matrix, this usually takes quite some time \n')
invcorrx = inv(corrx);
vif = diag(invcorrx);