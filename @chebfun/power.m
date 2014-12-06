function g = power(f, b, pref)
%.^   CHEBFUN power.
%   F.^G returns a CHEBFUN F to the scalar power G, a scalar F to the CHEBFUN
%   power G, or a CHEBFUN F to the CHEBFUN power G. F and or G may be complex.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT, COMPOSE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Initialise an empty CHEBFUN g:
g = chebfun();

% Empty case.
if ( isempty(f) || isempty(b) )
    return
end

if ( nargin < 3 )
    pref = chebfunpref();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHEBFUN .^ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ( isa(f, 'chebfun') && isa(b, 'chebfun') )
    
    % Check the number of columns match:
    if ( numColumns(f) ~= numColumns(b) )
        error('CHEBFUN:CHEBFUN:power:quasi', ...
            'Chebfun quasimatrix dimensions must agree.')
    end
    
    if ( numel(f) == 1 && numel(b) == 1 )
        % Array-valued CHEBFUN case:
        g = columnPower(f, b, pref);
    else
        % QUASIMATRIX case:

        % Convert to a cell array:
        f = mat2cell(f);
        b = mat2cell(b);
        % Loop over the columns:
        for k = numel(f):-1:1
            g(k) = columnPower(f{k}, b{k}, pref);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHEBFUN .^ constant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif ( isa(f, 'chebfun') )
    
    numColsF = numColumns(f);
    numelF = numel(f);
    numelB = numel(b);
    
    % Different cases:
    if ( numColsF == 1 )
        % e.g., x.^[1 2 3]
        for k = numelB:-1:1
            g(k) = columnPower(f, b(k), pref);
        end
        if ( all((b > 0) & (round(b) == b)) )
            g = quasi2cheb(g);
        end
    elseif ( numelB == 1 )
        % e.g., [x sin(x)].^2
        if ( numel(f) == 1 )
            % If F is array-valued, we might kick back a quasimatrix (to handle
            % singularites). Hence, we must assign to g, rather than g(1).
            g = columnPower(f, b, pref);
        else
            for k = numelF:-1:1
                g(k) = columnPower(f(k), b, pref);
            end
        end
    elseif ( numelB == numColsF )
        % e.g., [x sin(x) exp(x)].^[1 2 3]
        f = mat2cell(f);
        for k = numColsF:-1:1
            g(k) = columnPower(f{k}, b(k), pref);
        end
    else
        error('CHEBFUN:CHEBFUN:power:dim', ...
            'Chebfun quasimatrix dimensions must agree.');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% constant .^ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    
    numColsB = numColumns(b);
    numelF = numel(f);
    numelB = numel(b);
    
    % Different cases:
    if ( numColsB == 1 )
        % e.g., [1 2 3].^x
        for k = numelF:-1:1
            g(k) = columnPower(f(k), b, pref);
        end
        g = quasi2cheb(g);
    elseif ( numelF == 1 )
        % e.g., 1.^[x sin(x)] 
        if ( numelB == 1 )
            g = columnPower(f, b, pref);
        else
            b = mat2cell(b);
            for k = numelB:-1:1
                g(k) = columnPower(f, b{k}, pref);
            end
        end
    elseif ( numelF == numColsB )
        % e.g., [1 2].^[x sin(x)]
        b = mat2cell(b);
        for k = numColsB:-1:1
            g(k) = columnPower(f(k), b{k}, pref);
        end
        if ( numelB == 1 )
            g  = quasi2cheb(g);
        end
    else
        error('CHEBFUN:CHEBFUN:power:dim', ...
            'Chebfun quasimatrix dimensions must agree.');
    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = columnPower(f, b, pref)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHEBFUN .^ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if ( isa(f, 'chebfun') && isa(b, 'chebfun') ) 
    
    % Handle the periodic case:
    if ( isPeriodicTech(f) && ~isPeriodicTech(b) )
        f = chebfun(f);
    end
    
    % Call COMPOSE(): (Note, COMPOSE() checks that the domains match)
    g = compose(f, @power, b, pref);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHEBFUN .^ constant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( isa(f, 'chebfun') )
    
    if ( b == 0 )                       % Trivial case
        % Constant CHEBFUN:
        vals = ones(1, numColumns(f));
        g = chebfun(vals, f.domain([1,end]));
        
    elseif ( b == 1 )                   % Identity
        g = f;
        
    elseif ( b == 2 )                   % Square
        % Call TIMES():
        g = f.*f;
       
    elseif ( ( b > 0 ) && ( round(b) == b ) )   % Positive integer
        
        % If SINGFUN is involved:
        if ( issing(f) )
            % If singfun is involved, treat each piece individually:
            g = f;
            numFuns = numel(f.funs);
            for k = 1:numFuns
                g.funs{k} = power(f.funs{k}, b);
            end
            g.pointValues = g.pointValues.^b;
        else
            % Result will be smooth. Call COMPOSE():
            g = compose(f, @(x) power(x, b));
        end
        
    else                                % General case (SQRT is included)
        
        if ( numColumns(f) > 1 )
            % TODO: SINGFUN cannot handle array-valued CHEBFUNs, so we must
            % convert to a quasimatrix.
            f = cheb2cell(f);
            for k = numel(f):-1:1
                f{k} = columnPower(f{k}, b, pref);
            end
            g = quasimatrix(f);
            return
        end
        
        % Add breaks at the appropriate roots of f:
        if ( isreal(f) || round(b) == b )
            f = addBreaksAtRoots(f);
        else
            % Add breaks at the roots of the imaginary part of F to account for
            % the discontinuity in POWER along the negative real semi-axis due 
            % to the branch cut.
            r = getRootsForBreaks(imag(f));
            f = addBreaks(f, r);
            
        end
        % Loop over each piece individually:
        numFuns = numel(f.funs);
        g = f;
        for k = 1:numFuns
            g.funs{k} = power(f.funs{k}, b);
        end
        % Update the pointValues:
        g.pointValues = g.pointValues.^b;
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% constant .^ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    % Call COMPOSE():
    g = compose(b, @(x) power(f, x), pref);

end

end
