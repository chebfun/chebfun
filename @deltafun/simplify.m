function f = simplify(f)
%SIMPLIFY   Removes trivial rows and columns of impluses of a DELTAFUN F.
%   SIMPLIFY(F) removes trivial rows and columns from the magnitude matrix of the 
%   DELTAFUN F based on the tolerance. 
%
% See also SUM, CUMSUM.


deltaLoc = f.delta.location;
deltaMag = f.delta.magnitude;
% Merge columns if location of deltafunction are almost equal:
[deltaMag, deltaLoc] = mergeColumns(deltaMag, deltaLoc);

% Remove trivial columns:
m = length(deltaLoc);
j = 1;
for k = 1:m
    if( max(abs(deltaMag(:, j))) < deltafun.pref.deltafun.deltaTol )
        deltaLoc(j) = [];
        deltaMag(:, j) = [];
    else
        j = j + 1;
    end
end

% Remove trivial rows:
while( max(abs(deltaMag(end,:))) < deltafun.pref.deltafun.deltaTol )
    deltaMag = deltaMag(1:end-1,:);
end

% If any of these is empty, make everything empty explicitly. This is to avoid 
% annoying cases of m x 0 empty matrix or 0 x n empty matrix.
if ( isempty(deltaLoc) || isempty(deltaMag) )
    deltaLoc = [];
    deltaMag = [];
end

% Assign back:
f.delta.location = deltaLoc;
f.delta.magnitude = deltaMag;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Merge Columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, v] = mergeColumns(A, v)
% Merge columns of A if two locations given in v are almost equal.

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