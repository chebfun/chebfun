function h = plus(f, g)
% + PLUS of two BALLFUNV. 
%   F + G if F and G are BALLFUNV does componentwise addition. 
%
%   F + G if F is a double and G is a BALLFUNV does componentwise addition. 
% 
%   F + G if F is a BALLFUNV and G is a double does componentwise addition.
% 
%   PLUS(F,G) is called for the syntax F + G. 
% 
% See also MINUS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fIsBallfunv = isa(f, 'ballfunv');
gIsBallfunv = isa(g, 'ballfunv');

% Empty check
if isempty( f )
    h = g;
    if ~gIsBallfunv
        h = ballfunv(h);
    end
    return
end

% Empty check
if isempty( g )
    h = f;
    if ~fIsBallfunv
        h = ballfunv(h);
    end
    return
end

if ( fIsBallfunv && gIsBallfunv )
    F = f.comp;
    G = g.comp;
    h = ballfunv(F{1}+G{1},F{2}+G{2},F{3}+G{3});
elseif ( fIsBallfunv && isnumeric(g) )
    F = f.comp;
    h = ballfunv(F{1}+g,F{2}+g,F{3}+g);
elseif ( isnumeric(f) && gIsBallfunv)
    G = g.comp;
    h = ballfunv(G{1}+f,G{2}+f,G{3}+f);
else
    error('BALLFUNV:plus:unknown', ...
          ['Undefined function ''plus'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end

end
