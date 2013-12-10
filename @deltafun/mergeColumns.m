function [A, v] = mergeColumns(A, v)
%MERGECOLUMNS Merge columns of A if two locations given in V are almost equal.
%
%

% Get the tolerance:
tol = deltafun.pref.deltafun.proximityTol;

m = size(A, 2);
if ( length(v) ~= m || size(v, 1) ~= 1 )
    error( 'CHEBFUN:DELTAFUN:mergeColumns', 'No. of columns of A should be equal to the length of the vector va' );
end

% Make sure the input is sorted:
[v, idx] = sort(v);
A = A(:, idx);

j = 2;
for k = 2:m
    % Duplication scenarios
    p = ( v(j) == 0 && v(j-1) == 0 );
    vcMax = max(abs([v(j), v(j-1)]));
    p = p | ( abs((v(j)-v(j-1)))/vcMax < tol );
    % If there is a duplicate
    if ( p )
        % Merge the two columns
        A(:, j-1) = A(:, j-1) + A(:, j);
        % Remove the copied column of A and the location in va
        A(:, j) = [];
        v(j) = [];
    else
        j = j + 1;
    end
end

end