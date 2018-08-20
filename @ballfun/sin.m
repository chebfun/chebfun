function g = sin(f)
% SIN Sinus of a BALLFUN function
%   SIN(f) is the sinus of the BALLFUN function f
F = f.coeffs;
G = ballfun.vals2coeffs(sin(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
