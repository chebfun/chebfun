function y = polyfit(y, n)
%POLYFIT   Fit polynomial to a CLASSICFUN.
%   F = POLYFIT(Y, N) returns a CLASSICFUN F corresponding to the polynomial of
%   degree N that fits the CLASSICFUN Y in the least-squares sense.
%
%   Note CLASSICFUN/POLYFIT does not not support more than one output argument in
%   the way that MATLAB/POLYFIT does.
%
%   If y is a global polynomial of degree n then this code has an O(n (log n)^2)
%   complexity. If y is piecewise polynomial then it has an O(n^2) complexity.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Simply POLYFIT() the /ONEFUN:
y.onefun = polyfit(y.onefun, n);

end
