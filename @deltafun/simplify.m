function f = simplify(f)
%SIMPLIFY  Simplifis a DELTAFUN object F.
%
%   SIMPLIFY(F) removes trivial rows and columns from the magnitude matrix of the 
%   DELTAFUN F based on the tolerance and merges columns in the impulse matrix if
%   the location of delta functions is really close.
%
% See also SUM, CUMSUM.


deltaLoc = f.location;
deltaMag = f.impulses;
% Merge columns if location of deltafunction are almost equal:
[deltaMag, deltaLoc] = deltafun.mergeColumns(deltaMag, deltaLoc);

% Remove columns which are entriely below tolerance.
[deltaMag, deltaLoc] = deltafun.cleanColumns(deltaMag, deltaLoc);

% Remove ending rows of zeros.
deltaMag = deltafun.cleanRows(deltaMag);


% If any of these is empty, make everything empty explicitly. This is to avoid 
% annoying cases of m x 0 empty matrix or 0 x n empty matrix.
if ( isempty(deltaLoc) || isempty(deltaMag) )
    deltaLoc = [];
    deltaMag = [];
end

% Assign back:
f.location = deltaLoc;
f.impulses = deltaMag;

end