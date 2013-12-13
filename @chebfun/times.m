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
    if ( ~g.isTransposed )
        f = times(g, f);
    else
        f = times(g.', f).';
    end

elseif ( isempty(g) )          % CHEBFUN * []

    f = [];
    return

elseif ( isnumeric(g) )        % CHEBFUN * double

    % Loop over the FUNS:
    for k = 1:numel(f.funs)
        f.funs{k} = times(f.funs{k}, g);
    end

    % Multiply the pointValues:
    f.pointValues = f.pointValues .* g;

elseif ( ~isa(g, 'chebfun') )  % CHEBFUN * ???

    error('CHEBFUN:times:unknown', ...
          ['Undefined function ''times'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));

elseif ( isempty(f) )          % empty CHEBFUN * CHEBFUN

    % Nothing to do here. (Return an empty CHEBFUN as output.)
    return

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

    % Multiply the pointValues:
    f.pointValues = f.pointValues .* g.pointValues;

end

% Set small breakpoint values to zero:
f = thresholdBreakpointValues(f);

end
