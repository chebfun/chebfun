function f = plus(f, g)
%+	  Plus.
%   F + G adds chebfuns F and G, or a scalar to a chebfun if either F or G is a
%   scalar.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(f, 'chebfun') )   % ??? * chebfun

    % Ensure chebfun is the first input:
    f = times(g, f);
    return

elseif ( isempty(g) )       % chebfun + []

    f = [];
    return

elseif ( isnumeric(g) )     % chebfun + double

    % Add g to the funs:
    for k = 1:numel(f.funs)
        f.funs{k} = f.funs{k} + g;
    end
    
    % Add g to the impulses:
    g = repmat(g, size(f.impulses, 1), size(f.impulses, 2)./length(g));
    f.impulses = f.impulses + g;
    
    return

elseif ( ~isa(g, 'chebfun') ) % chebfun * ???

    error('CHEBFUN:times:unknown', ...
        'Undefined function ''times'' for input arguments of type %s and %s.', ...
        class(f), class(g));

elseif ( isempty(f) )       % empty chebfun + chebfun

    return

end

% Only chebfun + chebfun remains.

% Check to see if one of the chebfuns is transposed:
if ( xor(f.isTransposed, g.isTransposed) )
    error('CHEBFUN:plus:matdim', ...
        'Matrix dimensions must agree. (One input is transposed).');
end

% Overlap the chebfun objects:
[f, g] = overlap(f, g);

% Add the impulses:
f.impulses = f.impulses + g.impulses;

% Add the funs:
for k = 1:numel(f.funs)
    f.funs{k} = f.funs{k} + g.funs{k};
end

end