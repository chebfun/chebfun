function h = dot(f, g)
%DOT  Dot product of two BALLFUNV
%   DOT(f, g) is the dot of the BALLFUNV f and g

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;
G = g.comp;
h = F{1}.*G{1}+F{2}.*G{2}+F{3}.*G{3};
end
