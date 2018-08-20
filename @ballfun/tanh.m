function g = tanh(f)
% TANH Hyperbolic tangent of a BALLFUN function
%   TANH(f) is the hyperbolic tangent of the BALLFUN function f
F = f.coeffs;
G = ballfun.vals2coeffs(tanh(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
