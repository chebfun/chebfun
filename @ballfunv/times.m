function h = times(f, g)
%.*   Times of two BALLFUNV objects. 
%   F.*G if F is a BALLFUNV and G is double returns the BALLFUNV after
%   componentwise multiplication.
%
%   F.*G if F is a BALLFUNV and G is a BALLFUNV returns the BALLFUNV
%   after multiplication of F by each component of G.
% 
% See also MTIMES.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f ) || isempty(g)
    h = ballfunv();
    return
end

F = f.comp;
G = g.comp;
h = [F{1}.*G{1} ; F{2}.*G{2} ; F{3}.*G{3}];
end
