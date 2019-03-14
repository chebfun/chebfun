function h = dot(f, g)
%DOT   Vector dot product.
%   DOT(F, G) returns the dot product of the BALLFUNV F and G. 
% 
% See also CROSS. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty(f) || isempty(g)
    h = ballfun();
    return
end

F = f.comp;
G = g.comp;
h = F{1}.*G{1}+F{2}.*G{2}+F{3}.*G{3};
end
