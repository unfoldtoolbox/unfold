function B = Bernstein(x, t, j, k, alpha, leftflag)
% B = Bernstein(x, t, j, k)
% B = Bernstein(x, t, j, k, alpha)
% Compute Bernstein polynomial basis using de Casteljau's algorithm
%
% INPUTS:
%   x: vectors/array, point coordinates at which the function is to be
%      evaluated
%   t: vector, knots points, must be ascending sorted
%   j: vector, vector of spatial index, must be in [1:length(t)-k]
%      if it's empty all the basis functions are computed.
%   k-1: "order" of the spline (k is scalar)
%      k==1 -> piecewise constant
%      k==2 -> linear
%      k==4 -> cubic
%   alpha: vectors of size length(t)-k, optional coefficients of the basis
% OUTPUTS:
%   B: (m x n) where m is length(x), n is length(j)              
%       Each column of B is the basis function B_j,k
%   If t is not monotonically non-decreasing, B will be an empty array
%
%   Note: B_j,k has support on [t(j),t(j+k)[
%
%  Call with TRUE for 6th argument: B = Bernstein(..., TRUE) to compute
%  the "Left B-spline", which has support on ]t(j),t(j+k)].
%  The left B-spline is needed e.g. to evaluate recursion integral between
%  two B-splines.
%
% Author: Bruno Luong
%   31-May-2010: correct bug for LEFT Bspline is called with scalar
%   10-Jun-2010: Change in Bernstein.m to avoid NaN for end knots
% 
% Copyright (c) 2009, Bruno Luong
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%%
coefin = nargin>=5 && ~isempty(alpha);

if coefin
    if size(alpha,1)==1
        alpha = alpha.';
    end
end

cls = class(x);

%%
if isvector(x) %&& ~coefin
    szx = length(x);
else
    szx = size(x);
end
x = x(:);
m = size(x,1);

%%
% Special case, we ignore the Dirac distribution
if k<=0    
    if coefin
        B = zeros([szx size(alpha,2)],cls);
    else
        B = zeros([szx numel(j)],cls);
    end
    return
end

%%
% Max possible value of indice
maxj = length(t)-k;
if isempty(j) || coefin
    % all the index j
    j = 1:maxj;
    js = j;
else
    js  = sort(j(:));
end

% left and right bracket
jmin = js(1);
jmax = js(end);
% Check
if jmin<1 || jmax>maxj
    error('BERNSTEIN: j must be within [%d,%d]', 1, maxj);
end

%% k=1: step functions (piecwise constant)
B = zeros(m,jmax+k-jmin, cls);

if nargin>=6 && leftflag
    % Left B-spline
    tt = t(jmax+k:-1:jmin);
    if issorted(-tt)
        [trash col] = histc(-x,-tt);
        col = length(tt)-col; % Correct BUG, 31/05/2010
    else
        B = [];
        return
    end
else
    % Default, right B-spline
    tt = t(jmin:jmax+k);
    if issorted(tt)
        [trash col] = histc(x,tt); %#ok
    else
        B = [];
        return
    end
end

row = find(col>=1 & col<=size(B,2));
col = col(row);
B(row+(col-1)*m) = 1;


%%
for kk=2:k
    jvec = jmin:jmax+k-kk;
    dtvec = t(jvec+kk-1)-t(jvec);
    % recursion
    for c=1:length(jvec);
        jj = jvec(c);
        dt = dtvec(c);
        if dt~=0
            w1 = (x-t(jj)) / dt;
        else
            w1 = zeros(size(x),cls);
        end
        dt = (t(jj+kk)-t(jj+1));
        if dt~=0
            w2 = (t(jj+kk)-x) / dt;
        else
            w2 = zeros(size(x),cls);
        end
        ij = jj-jmin+1;
        B(:,ij) = w1.*B(:,ij) + w2.*B(:,ij+1);
    end
end

if length(j) ~= size(B,2)
    % Map to original vector j
    %[tf loc] = ismemberc(j, jmin:jmax); %#ok
    loc = ismembc2(j, jmin:jmax);
    B = B(:,loc);
end

%%
if coefin
    % Compute function from the coefficients
    B = multMat(B, alpha); % Bug fix 10-Jun-2010
    B = reshape(B, [szx size(alpha,2)]);
else
    % Basis
    B = reshape(B, [szx numel(j)]);
end

end % Bernstein

% Make the matrix product
function Balpha = multMat(B, alpha)
% Balpha = B*alpha; without returning NaN
col = any(B,1);
Balpha = B(:,col)*alpha(col,:);
% Balpha = B*alpha;
end
