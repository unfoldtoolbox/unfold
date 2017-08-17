function [ h] = imagesc( varargin)
%IMAGESC Summary of this function goes here
%   Detailed explanation goes here
warning('using the squeeze-slow-coverversion')
if length(varargin)==1 && length(size(varargin{1})) >2
    
    if sum(size(varargin{1}) ~= 1) == 2
        
        varargin{1} = squeeze(varargin{1});
    end
    
end
curdir=cd;
cd(fileparts(which('pie3')))

try
    hh = imagesc(varargin{:});
catch
end
cd(curdir)
% hh = builtin('imagesc',varargin{:});


if nargout > 0
    h = hh;
end

end

