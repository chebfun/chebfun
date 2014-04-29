function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum of the DELTAFUN F.
%   VALS = MINANDMAX(F) returns a 2-vector VALS = [MIN(F); MAX(F)] with the
%   global minimum and maximum of the DELTAFUN F.
%
%   [VALS, POS] = MINANDMAX(F) returns also the 2-vector POS where the minimum
%   and maximum of F occur.
%
% See also MIN, MAX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~anyDelta(f) )
    
    % The function has no delta functions:
    [vals, pos] = minandmax(f.smoothPart);
    
else
    
    % If there are delta functions, maxima and minima are undefined:
    vals = [NaN; NaN];
    pos  = [NaN; NaN];
end
  
end
