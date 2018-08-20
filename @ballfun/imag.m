function g = imag(f)
% IMAG Imaginary part of a BALLFUN function
%   IMAG(f) is the imaginary part of the BALLFUN function f

F = f.coeffs;
% Compute the imaginary part of the values and return the corresponding array of
% coefficients
G = ballfun.vals2coeffs(imag(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
