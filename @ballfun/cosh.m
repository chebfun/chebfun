function g = cosh(f)
% COSH Hyperbolic cosine of a BALLFUN function
%   COSH(f) is the hyperbolic cosine of the BALLFUN function f
F = f.coeffs;
G = ballfun.vals2coeffs(cosh(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
