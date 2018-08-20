function g = sinh(f)
% SINH Hyperbolic sinus of a BALLFUN function
%   SINH(f) is the hyperbolic sinus of the BALLFUN function f
F = f.coeffs;
G = ballfun.vals2coeffs(sinh(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
