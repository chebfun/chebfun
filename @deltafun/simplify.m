function f = simplify(f, pref)
%SIMPLIFY  Simplifys a DELTAFUN object F.
%   SIMPLIFY(F) removes trivial rows and columns from the magnitude matrix of
%   the DELTAFUN F based on the tolerance and merges columns in the impulse
%   matrix if the location of delta functions is really close.
%
% See also SUM, CUMSUM.

if ( nargin < 2 || isempty(pref) )
    pref = chebpref();
end

deltaLoc = f.deltaLoc;
deltaMag = f.deltaMag;

% Merge columns if location of deltafunction are almost equal:
[deltaMag, deltaLoc] = deltafun.mergeColumns(deltaMag, deltaLoc, pref);

% Remove columns which are entriely below tolerance:
[deltaMag, deltaLoc] = deltafun.cleanColumns(deltaMag, deltaLoc, pref);

% Remove ending rows of zeros:
deltaMag = deltafun.cleanRows(deltaMag, pref);

% If any of these is empty, return just the funPart:
if ( isempty(deltaLoc) || isempty(deltaMag) )
    deltaLoc = [];
    deltaMag = [];
    f = f.funPart;
    return
end

% Assign back:
f.deltaLoc = deltaLoc;
f.deltaMag = deltaMag;

end