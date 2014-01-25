function [g, jumpEnd] = cumsum(f, k, dim, shift)
%CUMSUM   Indefinite integral of a DELTAFUN.
%   CUMSUM(F) is the indefinite integral of the DELTAFUN F. If F has no delta
%   functions, then G is just the CUMSUM of the funPart of F. Any derivatives of
%   delta functions are integrated by once. This is done by shifting the rows of
%   the DELTAMAG matrix once upward. In case there are delta functions in F, a
%   cell array of FUNS is returned. If F has a delta function at the right end
%   point, it's magnitude is returned in JUMPEND, other wiase JUMPEND is set to
%   zero. The value in JUMPEND can be used by highe classes to CUMSUM 
%   conacatenated DELTAFUNS.
% See also SUM

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: Document and implemet k, dim and shift

% Trivial case:
if ( isempty(f) )
    % Return an empty deltafun in this case:
    g = deltafun;
    return;
end

deltaMag = f.deltaMag;
deltaLoc = f.location;

% Get tolerance:
pref = chebpref();
deltaTol = pref.deltaPrefs.deltaTol;

if ( isempty(deltaLoc) || isempty(deltaMag) )
    g = cumsum(f.funPart);
    jumpEnd = 0;
else
    % Clean up delta functions:
    f = simplify(f);
     
    % If f does not have a funPart, construct one compatible with delta
    % function locations:
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
                
    % Determine locations where jumps are to be introduced:
    idx = abs(deltaMag(1,:)) >= deltaTol;
    jumpLocs = deltaLoc(idx);
    
    % Determine the value of jumps;
    jumpVals = deltaMag(1, idx);
    
    % Integrate and clean the deltaMag matrix:
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
        % There are no jumps and hence no jump at the end:
        jumpEnd = 0;
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
        
        % The last jump:
        jumpEnd = jumpVals(end);
        
        % If the only jumps are at the left or the right end point, i.e. there
        % is no delta function in the interior, integrate and leave. There are
        % three possibilities of this sort:
        nJumps = length(jumpLocs);
        p1 = nJumps == 1 && deltaLeft;
        p2 = nJumps == 1 && deltaRight;
        p3 = nJumps == 2 && deltaLeft && deltaRight;
        if ( p1 || p2 || p3)
            s = cumsum(f.funPart) + jumpVals(1);
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
    end
end

end