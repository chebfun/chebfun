function g = restrict(f, s)
%RESTRICT   Restrict a DELTAFUN to a subinterval.
%   RESCTRICT(F, S) returns a DELTAFUN that is restricted to the subinterval
%   [S(1), S(2)] of the domain of F.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   an array of DELTAFUN objects, where the entries hold F restricted to each of
%   the subintervals defined by S.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    g = f;
    return
end

if ( isempty(f.funPart) )
    % If the funPart is empty, then there is no constraint on the values in s:
    a = -inf;
    b = inf;
else
    dom = f.funPart.domain;
    a = dom(1); 
    b = dom(2);
end

% Check if s is actually a subinterval:
if ( (s(1) < a) || (s(end) > b) || (any(diff(s) <= 0)) )
    error('DELTAFUN:restrict:badinterval', 'Not a valid interval.')
elseif ( (numel(s) == 2) && all(s == [a, b]) )
    % Nothing to do here!
    return
end

if ( ~isempty(f.funPart) )
    restrictedFuns = restrict(f.funPart, s);
else
    restrictedFuns = [];
end

if ( length(s) == 2 )
    % Only restricting to one subinterval, this will return a DELTAFUN:   
    g = deltafun();
    g.funPart = restrictedFuns;
    if ( ~isempty(f.location) )
        idx = (f.location >= s(1)) & (f.location <= s(2));
        g.location = f.location(idx);
        g.deltaMag = f.deltaMag(:, idx);
        if ( isempty(g.location) )
            g.location = [];
            g.deltaMag = [];
        end
    end
else
    % Restricting to multiple subintervals, this returns a cell-array of 
    % DELTAFUN objects.
    
    % Create a cell to be returned.
    g = cell(1, numel(s) - 1);
    
    % Create an empty DELTAFUN:
    emptyDeltaFun = deltafun();
    
    % Loop over each of the new subintervals, make a DELTAFUN and store in 
    % the cell returned:
    for k = 1:(numel(s) - 1)
        gk = emptyDeltaFun;
        
        if ( ~isempty(restrictedFuns) )
            gk.funPart = restrictedFuns{k};
        else
            gk.funPart = [];
        end
        
        if ( ~isempty(f.location) )
            idx = (f.location >= s(k)) & (f.location <= s(k+1));
            gk.location = f.location(idx);
            gk.deltaMag = f.deltaMag(:, idx);
        end
        
        if ( isempty(gk.location) )
            gk.location = [];
            gk.deltaMag = [];
        end
        g{k} = gk;
    end
end

end