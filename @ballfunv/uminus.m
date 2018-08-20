function g = uminus(f)
% UMINUS BALLFUNV unary minus
%   UMINUS(f) is negation of the BALLFUNV f
F = f.comp;
g = ballfunv(-F{1},-F{2},-F{3});
end
