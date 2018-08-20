function h = plus(f, g)
% PLUS Addition of two BALLFUN functions
%   PLUS(f, g) is addition of the BALLFUN functions f and g

fIsBallfun = isa(f, 'ballfun');
gIsBallfun = isa(g, 'ballfun');

if (fIsBallfun && gIsBallfun)
    X = f.coeffs+g.coeffs;
    h = ballfun(X);
elseif (fIsBallfun && isnumeric(g))
    S = size(f.coeffs);
    X = f.coeffs;
    % Add the constant g
    X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) = X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) + g;
    h = ballfun(X);
elseif (isnumeric(f) && gIsBallfun)
    S = size(g.coeffs);
    X = g.coeffs;
    % Add the constant f
    X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) = X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) + f;
    h = ballfun(X);
else
    error('BALLFUN:mtimes:unknown', ...
          ['Undefined function ''plus'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end
end
