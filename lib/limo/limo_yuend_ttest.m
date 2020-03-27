function [Ty,diff,se,CI,p,tcrit,df]=limo_yuend_ttest(a,b,percent,alpha)

% function [Ty,diff,se,CI,p,tcrit,df]=limo_yuend_ttest(a,b,percent,alpha)
%
% Computes Ty (Yuen's T statistic) to compare the trimmed means of two
% DEPENDENT groups (subjects with the same number of electrodes and frames).
% 
% INPUTS:
%
% a & b are 3D matrices Electrodes x Frames x Subjects. 
% percent must be a number between 0 & 100.
%
%   Default values:
%   percent = 20;
%   alpha = 0.05; 
%
% OUTPUTS:
%
% Ty   = Yuen T statistics. Ty is distributed approximately as Student's t 
%        with estimated degrees of freedom, df.
% diff = difference between trimmed means of a and b
% se   = standard error
% CI   = confidence interval around the difference
% p    = p value
%        The p-value (p) is the probability of obtaining a t-value whose absolute value
%        is greater than Ty if the null hypothesis (i.e., tma-tmb = mu) is true.
%        In other words, p is the p-value for a two-tailed test of H0: tma-tmb=0;
% tcrit = 1-alpha/2 quantile of the Student's t distribution with adjusted
% df    = dregrees of freedom
%
% See Wilcox (2005), Introduction to Robust Estimation and Hypothesis
% Testing (2nd Edition), page 188-191 for a description of the Yuen
% procedure for dependent groups.
%
% See also LIMO_TTEST LIMO_YUEN_TTEST 

% -----------------------------
%  Copyright (C) LIMO Team 2010

% original version: GAR, University of Glasgow, Dec 2007
% 3D, standalone version: GAR, University of Glasgow, June 2010
% GAR fixed bug for covariance / se

if nargin<4;alpha=.05;end
if nargin<3;percent=20;end

if isempty(a) || isempty(b) 
    error('limo_yuen_ttest:InvalidInput', 'data vectors cannot have length=0');
end

if (percent >= 100) || (percent < 0)
    error('limo_yuen_ttest:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

if percent == 50
    error('limo_yuend_ttest:InvalidPercent', 'PERCENT cannot be 50, use a method for medians instead.');
end

% number of trials
na=size(a,3);
nb=size(b,3);

if na ~= nb
    error('limo_yuend_ttest:InvalidInput', 'a and b must have the same size.');
else
    n = na;
end

g=floor((percent/100)*n); % number of items to winsorize and trim
h=n-2.*g; % effective sample size after trimming

% winsorise a
asort=sort(a,3);
loval=repmat(asort(:,:,g+1),[1 1 g+1]);
hival=repmat(asort(:,:,n-g),[1 1 g+1]);
wa=a;
for r=1:size(a,1)
    for c=1:size(a,2)
        wa(r,c,wa(r,c,:)<=loval(r,c))=loval(r,c);
        wa(r,c,wa(r,c,:)>=hival(r,c))=hival(r,c);
    end
end
wva=var(wa,0,3);

% winsorise b
bsort=sort(b,3);
loval=repmat(bsort(:,:,g+1),[1 1 g+1]);
hival=repmat(bsort(:,:,n-g),[1 1 g+1]);
wb=b;
for r=1:size(b,1)
    for c=1:size(b,2)
        wb(r,c,wb(r,c,:)<=loval(r,c))=loval(r,c);
        wb(r,c,wb(r,c,:)>=hival(r,c))=hival(r,c);
    end
end
wvb=var(wb,0,3);

% yuen's estimate of standard errors for a and b
da=(n-1).*wva;
db=(n-1).*wvb;

% covariance of winsorized samples
wcov=zeros(size(a,1),size(a,2));
for E=1:size(a,1)
    for F=1:size(a,2)
        tmp=cov(squeeze(wa(E,F,:)),squeeze(wb(E,F,:)));
        wcov(E,F)=tmp(1,2);
    end
end

dab=(n-1).*wcov;

% trimmed means
ma=mean(asort(:,:,(g+1):(n-g)),3);
mb=mean(bsort(:,:,(g+1):(n-g)),3);

diff=ma-mb;

df=h-1;
se=sqrt( (da+db-2.*dab)./(h.*(h-1)) );

Ty=diff./se;

p=2*(1-tcdf(abs(Ty),df)); % 2-tailed probability

tcrit=tinv(1-alpha./2,df); % 1-alpha./2 quantile of Student's distribution with df degrees of freedom
 
CI(:,:,1)=diff-tcrit.*se; 
CI(:,:,2)=diff+tcrit.*se;

end


