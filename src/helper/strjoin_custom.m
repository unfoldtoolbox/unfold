function [out] = strjoin_custom(in,varargin)
if nargin == 2
    string = varargin{1};
else
    string = ' ';
end
if iscell(in) && length(in)>1
    out = strjoin(in,string);
elseif iscell(in)
    out = in{1};
elseif ischar(in)
    out = in;
else
    error('unknown input')
end
