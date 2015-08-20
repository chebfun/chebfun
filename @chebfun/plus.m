function f = plus(f, g)
%+   CHEBFUN plus.
%   F + G adds CHEBFUNs F and G, or a scalar to a CHEBFUN if either F or G is a
%   scalar.
%
%   H = PLUS(F, G) is called for the syntax 'F + G'.
%
%   The dimensions of F and G must be compatible. Note that scalar expansion is
%   _not_ supported if both F and G are CHEBFUN objects.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(f, 'chebfun') )   % ??? + CHEBFUN
    
    % Ensure CHEBFUN is the first input:
    f = plus(g, f);
    return
    
elseif ( isempty(g) )       % CHEBFUN + []
    
    f = [];
    return
    
elseif ( isnumeric(g) )     % CHEBFUN + double
    
    if ( numel(f) == 1 )
        % Array-valued case:
        
        % Transpose g if f is transposed:
        if ( f(1).isTransposed )
            g = g.';
        end
        
        % Add g to the FUNs:
        for k = 1:numel(f.funs)
            f.funs{k} = f.funs{k} + g;
        end
        % Add g to the pointValues:
        if ( (size(f.pointValues, 2) == 1) &&  (numColumns(f) > 1) )
            f.pointValues = repmat(f.pointValues, 1, size(g, 2)); % Allow expansion in f.
        end
        if ( size(g, 2) > 1 )
            g = repmat(g, length(f.domain), 1);             % Allow expansion in g.
        end
        f.pointValues = f.pointValues + g;
    else
        % Quasimatrix case:
        
        numCols = numel(f);
        % Promote g if required:
        if ( isscalar(g) )
            g = repmat(g, 1, numCols);
        elseif ( length(g) ~= numCols || min(size(g)) ~= 1 )
            error('CHEBFUN:CHEBFUN:plus:dims', 'Matrix dimensions must agree.');
        end
        % Transpose g if f is a row CHEBFUN:
        if ( f(1).isTransposed )
            g = g.';
        end
        % Loop over the columns:
        for k = 1:numCols
            f(k) = f(k) + g(k);
        end
        
    end
    
elseif ( ~isa(g, 'chebfun') ) % CHEBFUN + ???
    
    error('CHEBFUN:CHEBFUN:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class(f), class(g));
    
elseif ( isempty(f) )         % empty CHEBFUN + CHEBFUN
    
    % Nothing to do. (Return empty CHEBFUN as output).
    return
    
else                          % CHEBFUN + CHEBFUN
    
    % Check to see if one CHEBFUN is transposed:
    if ( xor(f(1).isTransposed, g(1).isTransposed) )
        error('CHEBFUN:CHEBFUN:plus:matdim', ...
            'Matrix dimensions must agree. (One input is transposed).');
    elseif ( numColumns(f) ~= numColumns(g) )
        error('CHEBFUN:CHEBFUN:plus:dims', 'Matrix dimensions must agree.')
    end
    
    if ( numel(f) == 1 && numel(g) == 1 )
        % CHEBFUN case:
        
        % If one of the two CHEBFUNs uses a PERIODICTECH reprensetation, 
        % cast it to a NONPERIODICTECH.
        if ( ~isPeriodicTech(f.funs{1}) && isPeriodicTech(g.funs{1}) )
            g = chebfun(g, g.domain, 'tech', get(f.funs{1}, 'tech'));
        elseif ( isPeriodicTech(f.funs{1}) && ~isPeriodicTech(g.funs{1}) )
            f = chebfun(f, f.domain, 'tech', get(g.funs{1}, 'tech'));
        end
        
        % Overlap the CHEBFUN objects:
        [f, g] = overlap(f, g);
        
        % Add the pointValues:
        f.pointValues = f.pointValues + g.pointValues;
        % Add the FUNs:
        for k = 1:numel(f.funs)
            f.funs{k} = f.funs{k} + g.funs{k};
        end

    else
        % QUASIMATRIX case:
        
        % Loop over the columns:
        f = cheb2cell(f);
        g = cheb2cell(g);
        for k = numel(f):-1:1
            h(k) = f{k} + g{k};
        end
        f = h;
    end

end

% Set small breakpoint values to zero:
f = thresholdBreakpointValues(f);

end
