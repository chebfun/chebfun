function h = times(f, g)
%.*   BALLFUNV multiplication.
%   F.*G multiplies F and G, where F and G may be BALLFUNV objects or scalars.
%
% See also MTIMES.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%   F.*G if F is a ballfunv and G is double returns the ballfunv after
%   componentwise multiplication.F = f.comp;
F = f.comp;
G = g.comp;
h = [F{1}.*G{1} ; F{2}.*G{2} ; F{3}.*G{3}];
end
