function g = real(f)
% REAL Real part of a BALLFUN function
%   REAL(f) is the real part of the BALLFUN function f
F = f.coeffs;
% Compute the real part of the values and return the corresponding array of
% coefficients
G = ballfun.vals2coeffs(real(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
