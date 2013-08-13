function f = times(f, g)
%.*   Chebfun multiplication.
%   F.*G multiplies F and G, where F and G may be chebfun objects or scalars.
%   If F and/or G is array-valued, the dimensions must match.
%
% See also MTIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Deal with the special cases:
if ( ~isa(f, 'chebfun') )   % ??? * chebfun

    % Ensure chebfun is the first input:
    f = times(g, f);
    return

elseif ( isempty(g) )       % chebfun * []

    f = [];
    return

elseif ( isnumeric(g) )     % chebfun * double

    % Loop over the funs:
    for k = 1:numel(f.funs)
        f.funs{k} = times(f.funs{k}, g);
    end

    % Multiply the impulses:
    f.impulses = f.impulses .* g;

    return

elseif ( ~isa(g, 'chebfun') ) % chebfun * ???

    error('CHEBFUN:times:unknown', ...
        'Undefined function ''times'' for input arguments of type %s and %s.', ...
        class(f), class(g));

elseif ( isempty(f) )       % empty chebfun * chebfun

    return

end

% Only chebfun .* chebfun remains.

% Check to see if one of the chebfuns is transposed:
if ( xor(f.isTransposed, g.isTransposed) )
    error('CHEBFUN:times:matdim', ...
        'Matrix dimensions must agree. (One input is transposed).');
end

% Overlap:
[f, g] = overlap(f, g);

% Loop over the funs:
for k = 1:numel(f.funs)
    f.funs{k} = times(f.funs{k}, g.funs{k});
end

% Multiply the impulses:
f.impulses = f.impulses .* g.impulses;

end