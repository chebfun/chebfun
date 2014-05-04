function [g, rJump] = cumsum(f, k, dim)
%CUMSUM   Indefinite integral of a DELTAFUN.
%   [G, RJUMP] = CUMSUM(F) is the indefinite integral of the DELTAFUN F. If F 
%   has no delta functions, then G is just the CUMSUM of the funPart of F
%   RJUMP is the value of the delta function a the at the right end point of f.
%   Any derivatives of delta functions are integrated by shifting the rows 
%   of the DELTAMAG matrix once upward. 
%
%   In case there are delta functions in F, a cell array of FUNS is returned. 
%   If F has a delta function at the right end point, it's magnitude is added
%   to the right end point value of the funPart of the last deltafun. The 
%   value in RVAL can be used by higher classes to CUMSUM concatenated
%   DELTAFUNS.
%
% See also SUM

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: Document and implement k, dim

% Trivial case:
if ( isempty(f) )
    % Return an empty deltafun in this case:
    g = deltafun;
    return;
end

deltaMag = f.deltaMag;
deltaLoc = f.deltaLoc;

% Get tolerance:
pref = chebfunpref();
deltaTol = pref.deltaPrefs.deltaTol;

if ( isempty(deltaLoc) || isempty(deltaMag) )
    g = cumsum(f.funPart);
    % The jump at the right end is zero in this case:
    rJump = 0;
else
    % Clean up delta functions:
    f = simplifyDeltas(f);
     
    % Determine locations where jumps are to be introduced:
    idx = abs(deltaMag(1,:)) >= deltaTol;
    jumpLocs = deltaLoc(idx);
    
    % Determine the value of jumps;
    jumpVals = deltaMag(1, idx);
    
    % Remove the first row of deltaMag matrix since it is being integrated:
    if ( size(deltaMag, 1) > 1 )
        deltaMag = deltaMag(2:end, :);
        [deltaMag, deltaLoc] = deltafun.cleanColumns(deltaMag, deltaLoc);
    else
        deltaMag = [];
        deltaLoc = [];
    end
    
    if ( isempty(jumpLocs) )
        % No delta functions, but check for higher order delta functions:
        if ( isempty(deltaMag) || isempty(deltaLoc) )
            g = cumsum(f.funPart);
        else
            g = deltafun(cumsum(f.funPart), deltaMag, deltaLoc);
        end
        % There are no delta functions, so the jump is zero:
        rJump = 0;
        return
    else
        % There are delta functions, which will introduce jumps. First take care
        % of the end points.
        dom = f.funPart.domain;
        
        % If there is no delta function at the left end point, introduce a jump of
        % size 0 at the beginning:
        if ( jumpLocs(1) > dom(1) )
            deltaLeft = 0;
            jumpVals = [0, jumpVals];
        else
            deltaLeft = 1;
        end
        
        % If there is no delta function at the right end point, introduce a jump of
        % size 0 at the end:
        if ( jumpLocs(end) < dom(2) )
            deltaRight = 0;
            jumpVals = [jumpVals, 0];
        else
            deltaRight = 1;
        end
                
        % If the only jumps are at the left or the right end point, i.e. there
        % is no delta function in the interior, integrate and leave. There are
        % three possibilities of this sort:
        nJumps = length(jumpLocs);
        p1 = nJumps == 1 && deltaLeft;
        p2 = nJumps == 1 && deltaRight;
        p3 = nJumps == 2 && deltaLeft && deltaRight;
        if ( p1 || p2 || p3)
            s = cumsum(f.funPart) + jumpVals(1);
            rJump = jumpVals(end);
            if ( isempty(deltaLoc) || isemtpy(deltaMag) )
                g = s;                
            else
                g = deltafun(s, deltaMag, deltaLoc);
            end            
            return
        end
        
        %% Cell Array case:
        
        % Calculate the cumulative jump vector, entries of this vector will be used
        % to off-set the cumsum by the correct value.
        cumJump = cumsum(jumpVals);
        
        % Add end points to existing break points:
        breakPts = sort(union(f.funPart.domain, jumpLocs));
        
        % Get a cell array of funs via RESTRICT:
        funParts = restrict(f.funPart, breakPts);
        
        % Initialize output:
        nfuns = numel(funParts);
        g = cell(1, nfuns);
        for k = 1:nfuns
            % Integrate the funPart and add the appropriate jump:
            fk = cumsum(funParts{k}) + cumJump(k);
            domk = fk.domain;
            
            % Check whether there are delta functions in between the end points of
            % the current fun:
            idx = (deltaLoc >= domk(1)) & (deltaLoc <= domk(2));
            lk = deltaLoc(idx);
            
            % Construct a DELTAFUN or a FUN:
            if ( any(idx) )
                % Split end point delta functions, if exactly at the end-point
                % and if this fun is not the last fun:
                if ( k <= nfuns-1 && lk(end) == domk(2) )
                    lastIdx = find(idx, 1, 'last');
                    deltaMag(:, lastIdx) = deltaMag(:, lastIdx)/2;
                end
                mk = deltaMag(:, idx);               
                g{k} = deltafun(fk, mk, lk);
            else
                g{k} = fk;
            end
        end
        % Record the jump at the right end:
        rJump = jumpVals(end);
    end
end

end
