function f = simplifyDeltas(f, pref)
%SIMPLIFYDELTAS   Simplifys delta functions of a DELTAFUN object.
%   F = SIMPLIFYDELTAS(F) removes trivial rows and columns from the magnitude 
%   matrix of the DELTAFUN F based on the tolerance and merges columns in the 
%   impulse matrix if the location of delta functions is really close.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 || isempty(pref) )
    pref = chebfunpref();
elseif ( ~isa(pref, 'chebfunpref') )
    % Something other than a chebfunpref structure is passed. Assume it is a
    % tolerance:
    deltaTol = pref;
    pref = chebfunpref();
    pref.deltaPrefs.deltaTol = deltaTol;
end

deltaLoc = f.deltaLoc;
deltaMag = f.deltaMag;

% Merge columns if locations of delta function are almost equal:
[deltaMag, deltaLoc] = deltafun.mergeColumns(deltaMag, deltaLoc, pref);

% Remove columns which are entirely below tolerance:
[deltaMag, deltaLoc] = deltafun.cleanColumns(deltaMag, deltaLoc, pref);

% Remove ending rows of zeros:
deltaMag = deltafun.cleanRows(deltaMag, pref);

% If any of these is empty, return just the funPart:
if ( isempty(deltaLoc) || isempty(deltaMag) )
    f = f.funPart;
    return
end

% Assign back:
f.deltaLoc = deltaLoc;
f.deltaMag = deltaMag;

end
