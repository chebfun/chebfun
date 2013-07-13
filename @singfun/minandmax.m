function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum of the SINGFUN F on [-1,1].
%   VALS = MINANDMAX(F) returns a 2-vector VALS = [MIN(F); MAX(F)] with the
%   global minimum and maximum of the SINGFUN F on [-1,1].  

%   [TODO]: What would the following mean for a SINGFUN:
%   If F is a array-valued CHEBTECH, VALS is a 2-by-N matrix, where N is the number of
%   columns of F.  VALS(1, K) is the global minimum of the Kth column of F on
%   [-1, 1], and VALS(2, K) is the global maximum of the same.
%
%   [VALS, POS] = MINANDMAX(F) returns also the 2-vector POS where the minimum
%   and maximum of F occur.
%
%   [TODO]: Does this make sense for a SINGFUN:
%    If F is complex-valued the absolute values are taken to determine extrema
%   but the resulting values correspond to those of the original function. That
%   is, VALS = FEVAL(F, POS) where [~, POS] = MINANDMAX(ABS(F)). (In fact,
%   MINANDMAX actually computes [~, POS] = MINANDMAX(ABS(F).^2), to avoid
%   introducing singularlities to the function).
%
% See also MIN, MAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( any(f.exponents) )   
    if ( f.exponents(1) )
        % if positive infinity
        if ( feval(f,-1) > 0 )
            maxF = inf;
            maxLoc = -1;
        end
        if ( feval(f,-1) < 0 )
            minF = -inf;
            minLoc = -1;
        end
    end
    if ( f.exponents(2) )
        f2 = feval(f,1);
    end
    maxF = max(f1, f2);
    minF = 
else
    [vals, pos] = minandmax(f.smoothPart);
end

