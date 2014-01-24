function [g, jumpVals, locations] = cumsum(f)
%CUMSUM   Indefinite integral of a DELTAFUN.
%   CUMSUM(F) is the indefinite integral of the DELTAFUN F. LOCATIONS is a
%   vector which indicates the locations of the delta functions only (not their
%   derivatives). JUMPVALS is the vector of (signed) magnitude of jumps, that 
%   should be introduced at these locations.
% See also SUM

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output:
jumpVals = [];
locations = [];
g = deltafun;

% Trivial case:
if ( isempty(f) )
    % Return an empty deltafun in this case:
    return;
end

deltaMag = f.deltaMag;
deltaLoc = f.location;

if ( isempty(deltaLoc) || isempty(deltaMag) )
    g = cumsum(f.funPart);    
else
    % Clean up delta functions:
    f = simplify(f);
    
    % If f does not have a funPart, construct one matching the delta
    % function domain:    
    if ( isempty(f.funPart) )
        a = deltaLoc(1);
        if ( length(deltaLoc) < 2 )
            a = a-1;
            b = a+2;
        else
            b = deltaLoc(end);
        end
        f.funPart = fun.constructor(0, [a, b]);
    end
        
    % Get tolerance:
    pref = chebpref();
    deltaTol = pref.deltaPrefs.deltaTol;
        
    % Determine locations where jumps are to be introduced:
    idx = abs(deltaMag(1,:)) >= deltaTol;    
    newBreaks = deltaLoc(idx);
    
    % Determine the value of jumps, and remove them from the delta function
    % matrix:
    jumpVals = deltaMag(1, idx);
    dom = f.funPart.domain;    
    % If there is no delta function at the left end point, introduce a jump of
    % size 0:
    if ( deltaLoc(1) > dom(1) )
        jumpVals = [0, jumpVals];
    end
    
    % If there is no delta function at the right end point, introduce a jump of
    % size 0:
    if ( deltaLoc(end) < dom(2) )
        jumpVals = [jumpVals, 0];
    end
    
    % Remove the delta functions that have been integrated and simplify the
    % deltafunction matrix:
    deltaMag(1, :) = 0*deltaMag(1,:);                
    [deltaMag, deltaLoc] = deltafun.cleanColumns(deltaMag, deltaLoc);
    
    % Calculate the cumulative jump vector, entries of this vector will be used
    % to off-set the cumsum by the correct value.
    cumJump = cumsum(jumpVals);
    
    % Add end points to existing break points:
    breakPts = sort(union(f.funPart.domain, newBreaks));
    
    % Get a cell array of funs:
    funParts = restrict(f.funPart, breakPts);
    
    % Initialize output:
    nfuns = numel(funParts);
    g = cell(1, nfuns);
    for k = 1:numel(funParts)
        % Integrate the funPart and add the appropriate jump:
        fk = cumsum(funParts{k}) + cumJump(k);
        domaink = fk.domain;
        
        % Check whether there are delta functions in between the end points of
        % the current fun:
        idx = (deltaLoc >= domaink(1)) & (deltaLoc <= domaink(2));
        lk = deltaLoc(idx);
        
        % Construct a DELTAFUN or a FUN:
        if ( isempty(idx) )
            g{k} = fk;
        else
            mk = deltaMag(:, idx);    
            g{k} = deltafun(fk, mk, lk);
        end        
    end    
end
