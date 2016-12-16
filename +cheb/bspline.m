function B = bspline(m)
%BSPLINE   B-splines as chebfuns.
%   BSPLINE(M) returns the first M B-splines, starting from a constant value 1
%   on [-1/2,1/2]. Each sucessive B-spline is smoother by one derivative and
%   supported on an interval larger by 1/2 in each direction.
%
%   The returned result is a chebmatrix.
%
%   Example:
%      B = cheb.bspline(4); hold on
%      for k=1:4, plot(B{k}), end
%      axis([-0.1 1.1 -2 2])

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

s = chebfun(1, [-.5 .5]);
B = chebmatrix(s);
for k = 1:m-1
    B(k+1) = conv(B{k}, s);
end

end
