function m = mean(f)
%MEAN    Average or mean value of a BALLFUN.
%   MEAN(F) is the mean of the BALLFUN F.
%
% See also SUM3.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

m = sum3(f)*3/(4*pi);
end
