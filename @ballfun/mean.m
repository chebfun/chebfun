function m = mean(f)
% MEAN Mean of a BALLFUN function on the ballfun
%   MEAN(f) is the mean of the BALLFUN function f on the ballfun

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

m = sum3(f)*3/(4*pi);
end
