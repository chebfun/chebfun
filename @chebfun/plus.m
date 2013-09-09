function f = plus(f, g)
%+   CHEBFUN plus.
%   F + G adds CHEBFUNs F and G, or a scalar to a CHEBFUN if either F or G is a
%   scalar.
%
%   H = PLUS(F, G) is called for the syntax 'F + G'.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

if ( ~isa(f, 'chebfun') )   % ??? + CHEBFUN

    % Ensure CHEBFUN is the first input:
    f = plus(g, f);
    return

elseif ( isempty(g) )       % CHEBFUN + []

    f = [];

elseif ( isnumeric(g) )     % CHEBFUN + double
    
    % Transpose g if f is transposed:
    if ( f.isTransposed )
        g = g.';
    end

    % Add g to the FUNs:
    for k = 1:numel(f.funs)
        f.funs{k} = f.funs{k} + g;
    end
    
    % Add g to the impulses:
    if ( (size(f.impulses, 2) == 1) &&  (min(size(f)) > 1) )
        f.impulses = repmat(f.impulses, 1, size(g, 2)); % Allow expansion in f.
    end
    if ( size(g, 2) > 1 )
        g = repmat(g, length(f.domain), 1);             % Allow expansion in g.
    end
    f.impulses(:,:,1) = f.impulses(:,:,1) + g;

elseif ( ~isa(g, 'chebfun') ) % CHEBFUN + ???

    error('CHEBFUN:plus:unknown', ...
          ['Undefined function ''plus'' for input arguments of type %s ' ...
           'and %s.'], class(f), class(g));

elseif ( isempty(f) )         % empty CHEBFUN + CHEBFUN

    % Nothing to do. (Return empty CHEBFUN as output).

else                          % CHEBFUN + CHEBFUN

    % Check to see if one CHEBFUN is transposed:
    if ( xor(f.isTransposed, g.isTransposed) )
        error('CHEBFUN:plus:matdim', ...
            'Matrix dimensions must agree. (One input is transposed).');
    end

    % Overlap the CHEBFUN objects:
    [f, g] = overlap(f, g);

    % Add the impulses:
    f.impulses = f.impulses + g.impulses;

    % Add the FUNs:
    for k = 1:numel(f.funs)
        f.funs{k} = f.funs{k} + g.funs{k};
    end

end

end
