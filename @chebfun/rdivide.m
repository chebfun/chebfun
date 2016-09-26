function h = rdivide(f, g, pref)
%./   Pointwise CHEBFUN right divide.
%   F./G returns a CHEBFUN that represents the function F(x)/G(x).
%   If F and G are array-valued column (row) CHEBFUNs, they must have the same
%   number of columns (rows).
%
% See also MRDIVIDE, TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
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

if ( nargin < 3 )
    pref = chebfunpref();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHEBFUN ./ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ( isa(f,'chebfun') && isa(g, 'chebfun') )

    % Check the number of columns match:
    if ( numColumns(f) ~= numColumns(g) )
        error('CHEBFUN:CHEBFUN:rdivide:quasi', ...
            'Chebfun quasimatrix dimensions must agree.')
    end
    
    if ( numel(f) == 1 && numel(g) == 1 )
        % CHEBFUN case :
        
        % If one of the two CHEBFUNs uses a PERIODICTECH reprensetation, 
        % cast it to a NONPERIODICTECH.
        if ( ~isPeriodicTech(f.funs{1}) && isPeriodicTech(g.funs{1}) )
            g = chebfun(g, g.domain, 'tech', get(f.funs{1}, 'tech'));
        elseif ( isPeriodicTech(f.funs{1}) && ~isPeriodicTech(g.funs{1}) )
            f = chebfun(f, f.domain, 'tech', get(g.funs{1}, 'tech'));
        end
        
        % Array-valued CHEBFUN case:
        h = columnRdivide(f, g, pref);
    else
        % QUASIMATRIX case:
        
        % Convert to a cell array:
        f = mat2cell(f);
        g = mat2cell(g);
        % Loop over the columns:
        for k = numel(f):-1:1
            h(k) = columnRdivide(f{k}, g{k}, pref);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHEBFUN ./ constant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif ( isa(f, 'chebfun') )

    numColsF = numColumns(f);
    numelF = numel(f);
    numelG = numel(g);
    
    % Different cases:
    if ( numColsF == 1 )
        % e.g., x./[1 2 3]
        h = columnRdivide(f, g, pref);
    elseif ( numelG == 1 )
        % e.g., [x sin(x)]./2
        for k = numelF:-1:1
            h(k) = columnRdivide(f(k), g, pref);
        end
    elseif ( numelG == numColsF )
        % e.g., [x sin(x) exp(x)]./[1 2 3]
        if ( numel(f) == 1 )
            h = columnRdivide(f, g, pref);
        else
            f = mat2cell(f);
            for k = numColsF:-1:1
                h(k) = columnRdivide(f{k}, g(k), pref);
            end
        end
    else
        error('CHEBFUN:CHEBFUN:rdivide:dim', ...
            'Chebfun quasimatrix dimensions must agree.');
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% constant ./ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
   
    numColsG = numColumns(g);
    numelF = numel(f);
    numelG = numel(g);
    
    % Different cases:
    if ( numColsG == 1 )
        % e.g., [1 2 3]./(2+x)
        for k = numelF:-1:1
            h(k) = columnRdivide(f(k), g, pref);
        end
        h = quasi2cheb(h);
    elseif ( numelF == 1 )
        % e.g., 2./(2+[x x abs(x)])
        g = mat2cell(g);
        for k = numColsG:-1:1
            h(k) = columnRdivide(f, g{k}, pref);
        end
        try h = quasi2cheb(h); catch, end
    elseif ( numelF == numColsG )
        % e.g., [1 2]./(2+[x sin(x)])
        g = mat2cell(g);
        for k = numColsG:-1:1
            h(k) = columnRdivide(f(k), g{k}, pref);
        end
        try h = quasi2cheb(h); catch, end
    else
        error('CHEBFUN:CHEBFUN:rdivide:dim', ...
            'Chebfun quasimatrix dimensions must agree.');
    end

end

end

function h = columnRdivide(f, g, pref)

% If g is numeric then call TIMES():
if ( isnumeric(g) )
    if ( any(g == 0) )
        % note g could be a vector, with some zero components
        % TODO:  Return identically Inf/NaN CHEBFUN instead?
        error('CHEBFUN:CHEBFUN:rdivide:columnRdivide:divisionByZero', ...
            'Division by zero.')
    end
    h = f.*(1./g);
    return
end

% Check for zero FUNs:
for k = 1:numel(g.funs)
    if ( iszero(g.funs{k}) )
        % TODO:  Return CHEBFUN with identically Inf/NaN FUN instead?
        error('CHEBFUN:CHEBFUN:rdivide:columnRdivide:divisionByZeroChebfun', ...
            'Division by CHEBFUN with identically zero FUN.');
    end
end

% Add breaks at the roots of g:
g = addBreaksAtRoots(g, pref);

if ( isa(f, 'chebfun') )
    
    % Check that the domains are the same:
    if ( ~domainCheck(f, g) )
        error('CHEBFUN:CHEBFUN:rdivide:columnRdivide:domain', ...
            'Inconsistent domains.');
    end
    
    % Check that the orientation is the same:
    if ( xor(f.isTransposed, g.isTransposed) )
        error('CHEBFUN:CHEBFUN:rdivide:columnRdivide:dim', ...
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
