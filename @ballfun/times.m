function h = times(f, g)
%.*   Pointwise multiplication for BALLFUN.
%   F.*G multiplies BALLFUN F and G. Alternatively F or G could be 
%   a double.
%
% See also MTIMES.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isa(g, 'ballfunv')
    h = times(g,f);
else
    h = compose(f, @times, g);
end
end
