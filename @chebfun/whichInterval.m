function indx = whichInterval(dom, x)
%WHICHINTERVAL   Determine which interval a point lies in.
%   INDX = WHICHINTERVAL(DOM, X) returns a matrix of size(X) whos j,k entry is a
%   positive integer denoting which subinterval of the domain DOM (which should
%   be a sorted real-valued vector) the real part of the point X(j,k) is
%   positioned. INDX(j,k) = +/-INF if X(j,k) < DOM(1) or X(j,k) > DOM(end),
%   respectively.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

numInts = numel(dom) - 1;
xReal = real(x);
indx = NaN(size(x));

% Points to the left of the domain:
indx(xReal < dom(1)) = -inf;

% Points within the domain:
for j = 1:numInts
    indx( ( xReal >= dom(j) ) & ( xReal < dom(j+1) ) ) = j;
end
indx(xReal == dom(end)) = numInts;

% Points to the right of the domain:
indx(xReal > dom(end)) = -inf;

end
