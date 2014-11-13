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

% Get the domain:
dom = domain(f);
a = dom(1); 
b = dom(2);

% If the first or the last break points is very close to the original domain,
% merge it:
pref = chebfunpref();
tol = pref.deltaPrefs.proximityTol;
if ( abs(s(1) - a) < tol )
    s(1) = a;
end
if ( abs(s(end) - b) < tol )
    s(end) = b;
end
% Check if s is actually a subinterval:
if ( (s(1) < a) || (s(end) > b) || (any(diff(s) <= 0)) )
    error('CHEBFUN:DELTAFUN:restrict:badInterval', 'Not a valid subinterval.')
elseif ( (numel(s) == 2) && all(s == [a, b]) )
    % Nothing to do here!
    g = f;
    return
end

% Restrict the funPart:
restrictedFunParts = restrict(f.funPart, s);
if ( ~iscell(restrictedFunParts) )
    restrictedFunParts = {restrictedFunParts};
end

% Create a cell to be returned.
g = cell(1, numel(s)-1);

% Loop over each of the new subintervals, make a DELTAFUN and store in a cell:
numFuns = numel(s) - 1;
for k = 1:numFuns
    funPart = restrictedFunParts{k};
          
    idx = (f.deltaLoc >= s(k)) & (f.deltaLoc <= s(k+1));
    deltaLoc = f.deltaLoc(idx);
    deltaMag = f.deltaMag(:,idx);
    
    % If there are delta functions at a break point, divide them
    % by two. Each adjacent fun will get half the contribution. The first 
    % break point of the first fun and the last break point of the last 
    % fun do not get divided by half.
    if ( ~isempty(deltaLoc) )       
        if ( deltaLoc(1) == s(k) )
            if ( k ~= 1 )
                deltaMag(:, 1) = deltaMag(:, 1)/2;
            end
        end
                
        if ( deltaLoc(end) == s(k+1) )
            if ( k ~= numFuns )
                deltaMag(:, end) = deltaMag(:, end)/2;
            end
        end
    end
    
    if ( isempty(deltaLoc) )
        deltaLoc = [];
        deltaMag = [];
    end

    % Construct the new deltafun:
    if ( isempty(deltaLoc) )
        g{k} = funPart;
    else
        data.deltaMag = deltaMag;        
        data.deltaLoc = deltaLoc;
        g{k} = deltafun(funPart, data);
    end
end

% Return a DELTAFUN or CLASSICFUN if only one subinterval is requested:
if ( numel(s) == 2 )
    g = g{1};
end


end
