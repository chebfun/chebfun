function f = fun2ballfun(g,S)
% FUN2BALLFUN Conversion of a FUNCTION_HANDLE to a BALLFUN function
%   FUN2BALLFUN(g, S) is the BALLFUN function 
%   f = g(r, lambda, theta)

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Convert a handle_function to a ballfun function
m = S(1); n = S(2); p = S(3);

% Build the grid of evaluation points
r = chebpts(m);
lam = pi*trigpts(n);
th = pi*trigpts(p);

% Evaluate function handle at tensor grid:
[rr, ll, tt] = ndgrid(r, lam, th);
F = feval(g,rr, ll, tt); 

% Test if the function is constant
if size(F) == 1
   F = F(1)*ones(S);
end

f = ballfun(ballfun.vals2coeffs(F));
end
