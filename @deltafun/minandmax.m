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
    [vals, pos] = minandmax(f.funPart);
    
else
    % If there are delta functions, maxima and minima are inf or -inf:
    f = simplifyDeltas(f);
    deltaFuns = f.deltaMag(1, :);
    deltaLoc = f.deltaLoc;
    if ( isempty(deltaFuns) )
        % This means higher order deltas, return NaNs
        vals = [NaN; NaN];
        % Return the location of the first occurance of a higher order delta:
        pos  = [deltaLoc(1); deltaLoc(1)];
    else
        % This means at least one delta function so maxima and minima are inf 
        % and -inf regardless of the sign of the delta function.
        vals = [-inf, inf];        
        posIdx = find(deltaFuns > 0);
        negIdx = find(deltaFuns < 0);
        
        if ( isempty(negIdx) )
            % If there is no delta function with negative magnitude, return
            % the location of the first delta function.
            locMin = deltaLoc(1);
        else
            % If there is a delta function with negative magnitude, return the
            % location of first such occurance.
            locMin = deltaLoc(negIdx);
            locMin = locMin(1);
        end
        
        if ( isempty(posIdx) )
            % If there is no delta function with positive magnitude, return
            % the location of the first delta function.            
            locMax = deltaLoc(1);
        else
            % If there is a delta function with positive magnitude, return the
            % location of first such occurance.            
            locMax = deltaLoc(posIdx);
            locMax = locMax(1);
        end
        pos  = [locMin; locMax];
    end
end
  
end
