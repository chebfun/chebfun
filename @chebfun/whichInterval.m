function indx = whichInterval(dom, x, direction)
%WHICHINTERVAL   Determine which interval a point lies in.
%   INDX = WHICHINTERVAL(DOM, X) returns a matrix of size(X) whos j,k entry is a
%   positive integer denoting which subinterval of the domain DOM (which should
%   be a sorted real-valued vector) the real part of the point X(j,k) is
%   positioned. INDX(j,k) = -/+INF if X(j,k) < DOM(1) or X(j,k) > DOM(end),
%   respectively.
%   
%   WHICHINTERVAL(DOM, X, DIRECTION) determines which interval those points in X
%   which lie on the breakpoints of DOM should be counted in. DIRECTION = -1
%   (the default) places such points in the interval to the left, whilst
%   DIRECTION = 1 places them to the right.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numInts = numel(dom) - 1;
xReal = real(x);
indx = NaN(size(x));

% Points to the left of the domain:
indx(xReal < dom(1)) = -inf;

% Deal with the points that lie on breakpoints:
if ( nargin > 2 && direction > 0 )
    mygt = @(x, y) x >= y;
    mylt = @(x, y) x < y;
else
    mygt = @(x, y) x > y;
    mylt = @(x, y) x <= y;
end

indx(xReal == dom(1)) = 1;
% Points within the domain:
for j = 1:numInts
    indx( mygt(xReal, dom(j)) & mylt(xReal, dom(j+1)) ) = j;
end
indx(xReal == dom(end)) = numInts;

% Points to the right of the domain:
indx(xReal > dom(end)) = inf;

end
