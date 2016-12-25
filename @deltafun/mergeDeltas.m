function [D, w] = mergeDeltas( A, u, B, v)
%MERGEDELTAS   Merges delta functions and their locations.
%   [D, w] = mergeDeltas(A, u, B, v) assumes that A and B are matrices
%   representing delta functions and their derivatives at locations indicated by
%   the vectors U and V respectively. The function then tries to see if there
%   are locations in U and V that are the same or are very close to each other
%   and merges the corresponding delta functions in D and their locations in the
%   vector W.
%
% See also SIMPLIFY, CLEANROWS, CLEANCOLUMNS

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the sizes of both matrices:
szA = size(A);
szB = size(B);

%%
% Trivial cases:
if ( any(szA == 0 ) && any(szB == 0 ) )
    D = [];
    w = [];
    return
elseif ( any(szA == 0 ) )
    D = B;
    w = v;
    return
elseif ( any(szB == 0 ) )
    D = A;
    w = u;
    return
end

%%
% Both matrices are non empty, get the difference of rows:
d = szA(1) - szB(1);
if ( d >= 0 )
    % A has more rows, expand B:
    B = [B ; zeros(d, szB(2))];
else
    % B has more rows, expand A:
    A = [A ; zeros(-d, szA(2))];
end

% Concatenate:
D = [ A, B ];
w = [ u, v ];

% If w has duplicate locations, merge them:
[D, w] = deltafun.mergeColumns(D, w);
% If a column in D is entirely zero, remove it:
[D, w] = deltafun.cleanColumns(D, w);
% If there are trivial trailing rows, remove them:
D = deltafun.cleanRows(D);

end
