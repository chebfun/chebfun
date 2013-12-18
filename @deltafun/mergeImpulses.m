function [D, w] = mergeImpulses( A, u, B, v)
%MERGEIMPULSES   merges impulses matrices and their locations.
%   A and B are matrices. u

% Get the sizes of both matrices:
szA = size(A);
szB = size(B);

%%
% Trivial cases:
if ( any(szA == 0 ) && any(szB == 0 ) )
    D = [];
    w = [];
    return
end

if ( any(szA == 0 ) )
    D = B;
    w = v;
    return
end

if ( any(szB == 0 ) )
    D = A;
    w = u;
    return
end

%%
% Both matrices are non empty, get the differnece of rows:
d = szA(1) - szB(1);
if ( d >= 0 )
    % A has more rows, expand B:
    B = [ B; zeros(d, szB(2)) ];
else
    % B has more rows, expand A:
    A = [ A; zeros(-d, szA(2)) ];
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