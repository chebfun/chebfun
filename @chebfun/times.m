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
    f = times(g, f);

elseif ( isempty(g) )          % CHEBFUN * []

    f = [];

elseif ( isnumeric(g) )        % CHEBFUN * double

    % Loop over the funs:
    for k = 1:numel(f.funs)
        f.funs{k} = times(f.funs{k}, g);
    end

    % Multiply the impulses:
    f.impulses = f.impulses .* g;

elseif ( ~isa(g, 'chebfun') )  % CHEBFUN * ???

    error('CHEBFUN:times:unknown', ...
          ['Undefined function ''times'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));

elseif ( isempty(f) )          % empty CHEBFUN * CHEBFUN

    % Nothing to do here. (Return an empty CHEBFUN as output.)
    
else                           % CHEBFUN .* CHEBFUN 

    % Check to see if one of the CHEBFUNs is transposed:
    if ( xor(f.isTransposed, g.isTransposed) )
        error('CHEBFUN:times:matdim', ...
            'Matrix dimensions must agree. (One input is transposed).');
    end

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

end

end
