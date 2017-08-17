function [i] = get_min(a,B)
% This function returns the index of the minimum distance of all elemnts in
% B from a.
% If a is a vector, the function recursively calls itself.
% Keep in mind that multiple values could be the closest, this function
% selects the first one.
% > get_min([1 3],[5 4 3 2 1 2 3 4])
% > ans = [5 3]
if length(a) ==1
    [~,i] = min(abs(a-B));
    
else
    i = arrayfun(@(x)get_min(x,B),a);
end

