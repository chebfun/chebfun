function f = times(f, g)
%.*   CHEBFUN multiplication.
%   F.*G multiplies F and G, where F and G may be CHEBFUN objects or scalars.
%   If F and/or G is array-valued, the dimensions must match.
%
% See also MTIMES.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with the special cases:
if ( ~isa(f, 'chebfun') )      % ??? * CHEBFUN

    % Ensure CHEBFUN is the first input:
    if ( ~g(1).isTransposed )
        f = times(g, f);
    else
        f = times(g.', f.').';
    end

elseif ( isempty(g) )          % CHEBFUN * []

    f = [];
    return

elseif ( isempty(f) )          % empty CHEBFUN * CHEBFUN

    % Nothing to do here. (Return an empty CHEBFUN as output.)
    return    

elseif ( isnumeric(g) )        % CHEBFUN * double

    if ( numel(f) == 1 )
        % Array-valued case:

        % Loop over the funs:
        for k = 1:numel(f.funs)
            f.funs{k} = times(f.funs{k}, g);
        end

        % Multiply the pointValues:
        if ( numel(g) > 1 )
            f.pointValues = bsxfun(@times, f.pointValues, g);
        else
            f.pointValues = f.pointValues .* g;
        end
    
    else
        % Quasimatrix case:

        numCols = numel(f);
        % Promote g if required:
        if ( ~isscalar(g) )
            error('CHEBFUN:CHEBFUN:times:dim', 'Matrix dimensions must agree.');
        end
        % Loop over the columns:
        for k = 1:numCols
            f(k) = f(k).*g;
        end
        
    end


elseif ( ~isa(g, 'chebfun') )  % CHEBFUN * ???

    error('CHEBFUN:CHEBFUN:times:unknown', ...
          ['Undefined function ''times'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));

else                           % CHEBFUN .* CHEBFUN 

    % Check to see if one of the CHEBFUNs is transposed:
    if ( xor(f(1).isTransposed, g(1).isTransposed) )
        error('CHEBFUN:CHEBFUN:times:matdim', ...
            'Matrix dimensions must agree. (One input is transposed).');
    end
    
    if ( numColumns(f) ~= numColumns(g) )
        error('CHEBFUN:CHEBFUN:times:matdim', 'Matrix dimensions must agree.');
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
        
        % Overlap:
        [f, g] = overlap(f, g);

        % Loop over the FUNs:
        for k = 1:numel(f.funs)
            f.funs{k} = times(f.funs{k}, g.funs{k});
        end
        f.pointValues = f.pointValues .* g.pointValues;
        
    else
        % QUASIMATRIX case:

        % Loop over the columns:
        f = cheb2quasi(f);
        g = cheb2quasi(g);
        for k = 1:numel(f)
            f(k) = f(k).*g(k);
        end

    end
end

% Set small breakpoint values to zero:
f = thresholdBreakpointValues(f);

end
