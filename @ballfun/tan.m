function g = tan(f)
% TAN Tangent of a BALLFUN function
%   TAN(f) is the tangent of the BALLFUN function f
F = f.coeffs;
G = ballfun.vals2coeffs(tan(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
