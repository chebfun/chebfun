function [A, B] = merge(A, B)
%MERGE   Merge information from two CHEBDSICRETIZATION objects.
%   [A, B] = MERGE(A, B) merges two CHEBDISCRETOAZTIONS A and B.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% ENsure compatibility:
if ( class(A) ~= class(B) )
    error('CHEBFUN:chebDiscretization:merge:incompt', ...
        'Incompatible discretizations.');
end

% Merge the domains:
dom = chebfun.mergeDomains(A.domain, B.domain);
A.domain = dom;
B.domain = dom;

% Merge the discretization lengths:
dim = mymax(A.dimension, B.dimension);
A.dimension = dim;
B.dimension = dim;

% Merge the dimension adjustments:
dimAdjust = mymax(A.dimAdjust, B.dimAdjust);
A.dimAdjust = dimAdjust;
B.dimAdjust = dimAdjust;

% Merge the projection orders:
projOrder = mymax(A.projOrder, B.projOrder);
A.projOrder = projOrder;
B.projOrder = projOrder;

end

function out = mymax(a, b)
%MYMAX   Take maximum, but avoid empties.
if ( isempty(a) )
    out = b;
elseif ( isempty(b) )
    out = a;
else
    out = max(a, b);
end

end
