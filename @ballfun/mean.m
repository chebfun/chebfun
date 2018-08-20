function m = mean(f)
% MEAN Mean of a BALLFUN function on the ballfun
%   MEAN(f) is the mean of the BALLFUN function f on the ballfun
m = sum3(f)*3/(4*pi);
end
