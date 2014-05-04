function [A, v] = mergeColumns(A, v, pref)
%MERGECOLUMNS Merges columns of A if two locations given in V are almost equal.
%   [A v] = MERGECOLUMNS(A, V) merges two columns if the corresponding entries
%   in V are equal or close to each other more than a certain tolerance.
%
% See also CLEANCOLUMNS, MERGEDELTAS, CLEANROWS.

% Get the tolerance:
if ( nargin < 3 || isempty(pref) )
    pref = chebfunpref();
end
tol = pref.deltaPrefs.proximityTol;

m = size(A, 2);
if ( length(v) ~= m || size(v, 1) > 1 )
    error('CHEBFUN:DELTAFUN:mergeColumns', ...
        'Number of columns of A should be equal to the length of the vector v.');
end

% Make sure the input is sorted:
[v, idx] = sort(v);
A = A(:,idx);

j = 2;
for k = 2:m
    % Duplication scenarios.
    
    % If two entries are equal to zero:
    p = ( v(j) == 0 && v(j-1) == 0 );
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
    else
        j = j + 1;
    end
end

end