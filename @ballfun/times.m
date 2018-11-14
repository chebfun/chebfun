function h = times(f, g)
%.*   BALLFUN multiplication.
%   F.*G multiplies F and G, where F and G may be BALLFUN objects or scalars.
%   If F and/or G is array-valued, the dimensions must match.
%
% See also MTIMES.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = compose(f, @times, g); 
end
