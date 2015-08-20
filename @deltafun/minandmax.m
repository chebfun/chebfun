function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum of the DELTAFUN F.
%   VALS = MINANDMAX(F) returns a 2-vector VALS = [MIN(F); MAX(F)] with the
%   global minimum and maximum of the DELTAFUN F.
%
%   [VALS, POS] = MINANDMAX(F) returns also the 2-vector POS where the minimum
%   and maximum of F occur.
%
% See also MIN, MAX.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for the empty case:
if ( isempty(f) )
    vals = [];
    pos = [];
    return
end

% Take the minandmax of the funPart:
[vals, pos] = minandmax(f.funPart);
    
f = simplifyDeltas(f);
% Deal with the case when there are delta functions:
% NOTE: Higher order delta functions have no effect on maxima or minima of the
% function.
if ( isa(f, 'deltafun') )
    deltaMag = f.deltaMag;       % All delta functions in f.
    deltaFuns = deltaMag(1, :);  % zeroth order delta functions.
    deltaLoc = f.deltaLoc;   
    
    posIdx = find(deltaFuns > 0 );
    negIdx = find(deltaFuns < 0 );    
    
    % If there is a positive delta function, maxima is Inf:
    if ( ~isempty(posIdx) )
        vals(2) = Inf;
        loc = deltaLoc(posIdx);
        pos(2) = loc(1);
    end
    
    % If there is a negative delta function, minima is Inf:
    if ( ~isempty(negIdx) )
        vals(1) = -Inf;
        loc = deltaLoc(negIdx);
        pos(1) = loc(1);
    end        
end
  
end
