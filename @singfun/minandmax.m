function [vals, pos] = minandmax(f)
%MINANDMAX   Global minimum and maximum of the SINGFUN F on [-1,1].
%   VALS = MINANDMAX(F) returns a 2-vector VALS = [MIN(F); MAX(F)] with the
%   global minimum and maximum of the SINGFUN F on [-1,1].
%
%   [VALS, POS] = MINANDMAX(F) returns also the 2-vector POS where the minimum
%   and maximum of F occur.
%
% See also MIN, MAX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

tol = chebfunpref().blowupPrefs.exponentTol;

if ( ~any(f.exponents) || all(abs(f.exponents) < tol) )
  
    % The function is actually smooth!
    [vals, pos] = minandmax(f.smoothPart);
    
else
    
    % Initialise:
    minF = [];
    maxF = [];
    minLoc = [];
    maxLoc = [];   
    
    % Look for blow up at the left:
    if ( f.exponents(1) < -tol ) % Singularity at the left end.
        fVal = feval(f, -1);
        if ( fVal == inf )
            maxF = inf;
            maxLoc = -1;
        elseif ( fVal == -inf )
            minF = -inf;
            minLoc = -1;
        else
            % NaNs may occur and then we can not conclude anything.
            error('CHEBFUN:SINGFUN:minandmax:boundedness', ...
            ['Function has a singularity but isn''t infinite at the left ' ...
             'endpoint']);
        end
    end
    
    % Look for blow up at the right:
    if ( f.exponents(2) < -tol ) % Singularity at the right end.
        fVal = feval(f, 1);
        if ( fVal == inf )
            maxF = inf;
            maxLoc = 1;
        elseif ( fVal == -inf )
            minF = -inf;
            minLoc = 1;
        else
            % NaNs may occur and then we can not conclude anything.
            error('CHEBFUN:SINGFUN:minandmax:boundedness', ...
            ['Function has a singularity but isn''t infinite at the right ' ...
             'endpoint']);
        end
    end
    
    % If min or max is empty, then we need to do more work:
    if ( isempty(minF) || isempty(maxF) )       
        % Find the roots of the derivative for local minima:
        fp = diff(f);
        r = roots(fp);
        % Append the end points and remove duplicates:
        r = unique([-1 ; r ; 1]);
        if ( isempty(maxF) )
            % Take the maximum of the local maxima:
            [maxF, maxIndex] = max(feval(f, r));
            maxLoc = r(maxIndex);
        end
        if ( isempty(minF) )
            % Take the minimum of the local minima:
            [minF, minIndex] = min(feval(f, r));
            minLoc = r(minIndex);
        end             
    end    
    vals = [ minF; maxF ];
    pos = [ minLoc; maxLoc ];
    
end
  
end
