function g = cos(f)
% COS Cosine of a BALLFUN function
%   COS(f) is the cosine of the BALLFUN function f
F = f.coeffs;
G = ballfun.vals2coeffs(cos(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
