function h = rdivide(f, g)
%./   Pointwise CHEBFUN right divide.
%   F./G returns a CHEBFUN that represents the function F(x)/G(x).
%   If F and G are array-valued column (row) CHEBFUNs, they must have the same
%   number of columns (rows).
%
% See also MRDIVIDE, TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case:
if ( isempty(f) || isempty(g) )
    h = chebfun(); 
    return
end

% Trivial zero numerator case:
if ( isnumeric(f) && ~any(f) )
    h = 0*g;
    return
end

% If g is numeric then call TIMES():
if ( isnumeric(g) )
    if ( g == 0 )
        % TODO:  Return identically Inf/NaN CHEBFUN instead?
        error('CHEBFUN:rdivide:DivisionByZero', 'Division by zero.')
    end
    h = f.*(1./g);  
    return
end

% Check for zero FUNs:
for k = 1:numel(g.funs)
    if ( iszero(g.funs{k}) )
        % TODO:  Return CHEBFUN with identically Inf/NaN FUN instead?
        error('CHEBFUN:rdivide:DivisionByZeroChebfun', ...
            'Division by CHEBFUN with identically zero FUN.');
    end
end

% Add breaks at the roots of g:
g = addBreaksAtRoots(g);

if ( isa(f, 'chebfun') )
    
    % Check that the domains are the same:
    if ( ~domainCheck(f, g) )
        error('CHEBFUN:rdivide:domain', 'Inconsistent domains.');
    end
    
    % Check that the orientation is the same:
    if ( xor(f.isTransposed, g.isTransposed) )
        error('CHEBFUN:rdivide:dim', ...
            'Matrix dimension do not agree (transposed)');
    end
    
    % Introduce matching breakpoints in f and g:
    [f, g] = overlap(f, g);

    % Copy g to h in preparation for output:
    h = g;

    % Loop over the FUNS:
    for k = 1:numel(g.funs)
        h.funs{k} = rdivide(f.funs{k}, g.funs{k});
    end
    
    % Divide the pointValues:
    h.pointValues = f.pointValues./g.pointValues;
    
else

    % Copy g to h in preparation for output:
    h = g;

    % Loop over the FUNS:
    for k = 1:numel(g.funs)
        h.funs{k} = rdivide(f, g.funs{k});
    end
    
    % Divide the pointValues:
    h.pointValues = f./g.pointValues;
    
end

end
