function f = times(f, g)
%.*   CHEBFUN multiplication.
%   F.*G multiplies F and G, where F and G may be CHEBFUN objects or scalars.
%   If F and/or G is array-valued, the dimensions must match.
%
% See also MTIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

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

elseif ( isnumeric(g) )        % CHEBFUN * double

    if ( numel(f) == 1 )
        % Array-valued case:

        % Loop over the funs:
        for k = 1:numel(f.funs)
            f.funs{k} = times(f.funs{k}, g);
        end

        % Multiply the impulses:
        if ( numel(g) > 1 )
            f.impulses = bsxfun(@times, f.impulses, g);
        else
            f.impulses = f.impulses .* g;
        end
    
    else
        % Quasimatrix case:

        numCols = numel(f);
        % Promote g if required:
        if ( ~isscalar(g) )
            error('CHEBFUN:times:dim', 'Matrix dimensions must agree.');
        end
        % Loop over the columns:
        for k = 1:numCols
            f(k) = f(k).*g;
        end

    end


elseif ( ~isa(g, 'chebfun') )  % CHEBFUN * ???

    error('CHEBFUN:times:unknown', ...
          ['Undefined function ''times'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));

elseif ( isempty(f) )          % empty CHEBFUN * CHEBFUN

    % Nothing to do here. (Return an empty CHEBFUN as output.)
    
else                           % CHEBFUN .* CHEBFUN 

    % Check to see if one of the CHEBFUNs is transposed:
    if ( xor(f(1).isTransposed, g(1).isTransposed) )
        error('CHEBFUN:times:matdim', ...
            'Matrix dimensions must agree. (One input is transposed).');
    end

    if ( numel(f) == 1 && numel(g) == 1 )
        % Array-valued CHEBFUN case:

        % Overlap:
        [f, g] = overlap(f, g);

        % Loop over the FUNs:
        for k = 1:numel(f.funs)
            f.funs{k} = times(f.funs{k}, g.funs{k});
        end

        % Multiply the impulses:
        % [TODO]:  This doesn't make sense for higher-order impulses:  you can't
        % multiply two Dirac deltas!  What to do?
        f.impulses = f.impulses .* g.impulses;

    else
        % QUASIMATRIX case:

        if ( numel(f) ~= numel(g) )
            error('CHEBFUN:plus:dims', 'Matrix dimensions must agree.');
        else
            % Loop over the columns:
            for k = 1:numel(f)
                f(k) = f(k).*g(k);
            end
        end

    end


end

end
