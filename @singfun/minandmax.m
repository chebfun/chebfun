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
tol = singfun.pref.singfun.eps;
if ( any(f.exponents) )
    minF = [];
    maxF = [];
    minLoc = [];
    maxLoc = [];   
    if ( f.exponents(1) < -tol )
        % If there is a singularity at the left end
        fVal = feval(f, -1);
        if ( fVal == inf )
            maxF = inf;
            maxLoc = -1;
        elseif ( fVal == -inf )
            minF = -inf;
            minLoc = -1;
        else
            error('CHEBFUN:SINGFUN:minandmax:boundedness', ...
        'function has a singularity but bounded at the left end point');
        end
    end
    
    if ( f.exponents(2) < -tol )
        % if there is a singularity at the right end
        fVal = feval(f, 1);
        if ( fVal == inf )
            maxF = inf;
            maxLoc = 1;
        elseif ( fVal == -inf )
            minF = -inf;
            minLoc = 1;
        else
            error('CHEBFUN:SINGFUN:minandmax:boundedness', ...
        'function has a singularity but bounded at the right end point');
        end
    end
    
    % if min or max is empty, then we need to do more work
    if ( isempty(minF) || isempty(maxF) )       
        % find the roots of the derivative for local minima
        r = roots(diff(f));
        % append the end points and remove duplicates
        r = unique([-1;r;1]);
        if ( isempty(maxF) )
            % take the maximum of the local maxima
            [maxF, maxIndex] = max(feval(f,r));            
            maxLoc = r(maxIndex);
        end
        if ( isempty(minF) )
            % take the minimum of the local minima
            [minF, minIndex] = min(feval(f,r));          
            minLoc = r(minIndex);
        end             
    end    
    vals = [ minF; maxF ];
    pos = [ minLoc; maxLoc ];
else
    % the function is actually smooth. 
    [vals, pos] = minandmax(f.smoothPart);
end
  
end


