function g = log(f)
% LOG Logarithm of a BALLFUN function
%   LOG(f) is the logarithm of the BALLFUN function f

% Return the logarithm of the ballfun function f
F = f.coeffs;
G = ballfun.vals2coeffs(log(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
