function g = uminus(f)
% UMINUS BALLFUN unary minus
%   UMINUS(f) is negation of the BALLFUN function f
X = f.coeffs;
g = ballfun(-X);

end
