function [A, v, I] = mergeColumns(A, v, pref)
%MERGECOLUMNS Merges columns of A if locations given in V are almost equal.
%   [A, v, I] = MERGECOLUMNS(A, V) merges two columns if the corresponding
%   entries in V are equal or close to each other more than a certain
%   tolerance. I contains the indices where duplication existed.
%
% See also CLEANCOLUMNS, MERGEDELTAS, CLEANROWS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get the tolerance:
if ( nargin < 3 || isempty(pref) )
    pref = chebfunpref();
end
tol = pref.deltaPrefs.proximityTol;

m = size(A, 2);
if ( length(v) ~= m || size(v, 1) > 1 )
    error('CHEBFUN:DELTAFUN:mergeColumns:mergeColumns', ...
        'Number of columns of A should be equal to the length of the vector v.');
end

% Make sure the input is sorted:
[v, idx] = sort(v);
A = A(:,idx);
I = [];

j = 2;
for k = 2:m
    % Duplication scenarios.
    
    % If two entries are equal to zero:
    p = ( v(j) == 0 && v(j-1) == 0 );
    % Or if one of them is zero and the other close to zero or if both are close
    % to zero:
    if ( all(abs(v(j-1:j)) < 10*eps ) )
        p = p | true;
    end
        
    % Or if they are very close:
    vcMax = max(abs([v(j), v(j-1)]));
    p = p | ( abs((v(j)-v(j-1)))/vcMax < tol );
    
    
    % If there is a duplicate:
    if ( p )
        % Merge the two columns:
        A(:, j-1) = A(:, j-1) + A(:, j);
        % Remove the copied column of A and the corresponding location in v:
        A(:, j) = [];
        v(j) = [];
        I = [I, k]; %#ok<AGROW>
    else
        j = j + 1;
    end
end

end
