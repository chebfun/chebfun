% Test file for trigtech/vals2coeffs.m

function pass = test_vals2coeffs(varargin)

% Set a tolerance (pref.chebfuneps doesn't matter)
tol = 100*eps;

%%
% Test that a single value is converted correctly
vals = sqrt(2);
c = trigtech.vals2coeffs(vals);
pass(1) = (vals == c);

%%
% Simple data (odd case)
% Exact values
f_exact = @(x) 1 + cos(2*pi*x);
x_pts = trigtech.trigpts(5);
vals = f_exact(x_pts);
% Exact coefficients
cTrue = [0.5 0 1 0 0.5].';

%%
% Test real branch
c = trigtech.vals2coeffs(vals);
pass(2) = norm(c - cTrue, inf) < tol;

%%
% Test imaginary branch
c = trigtech.vals2coeffs(1i*vals);
pass(3) = norm(c - 1i*cTrue, inf) < tol;

%%
% Test general branch
c = trigtech.vals2coeffs((1+1i)*vals);
pass(4) = norm(c - (1+1i)*cTrue, inf) < tol;

%%
% Test for array input
c = trigtech.vals2coeffs([vals, -vals]);
pass(5) = norm(c(:,1) - cTrue, inf) < tol && ...
          norm(c(:,2) + cTrue, inf) < tol;
      
%%
% Simple data (even case)
x_pts = trigtech.trigpts(6);
% Exact values
f_exact = @(x) 1 + sin(2*pi*x) + cos(3*pi*x);
vals = f_exact(x_pts);
% Exact coefficients
cTrue = [1 0.5i 0 1 0 -0.5i].';

%%
% Test real branch
c = trigtech.vals2coeffs(vals);
pass(6) = norm(c - cTrue, inf) < tol;

%%
% Test imaginary branch
c = trigtech.vals2coeffs(1i*vals);
pass(7) = norm(c - 1i*cTrue, inf) < tol;

%%
% Test general branch
c = trigtech.vals2coeffs((1+1i)*vals);
pass(8) = norm(c - (1+1i)*cTrue, inf) < tol;

%%
% Test for array input
c = trigtech.vals2coeffs([vals, -vals]);
pass(9) = norm(c(:,1) - cTrue, inf) < tol && ...
          norm(c(:,2) + cTrue, inf) < tol;
      
%%
% Test for symmetry
x = trigpts(123);
vals = [cos(pi*x),sin(pi*x),cos(pi*x)+sin(pi*x)];
vals(1,2) = 0;
c = trigtech.vals2coeffs(vals);
pass(10) = norm(c(:,1) - flipud(c(:,1)), inf) == 0 && ...
           norm(c(:,2) + flipud(c(:,2)), inf) == 0 && ...
           norm(c(:,3) - flipud(c(:,3)), inf) > 0 && ...
           norm(c(:,3) + flipud(c(:,3)), inf) > 0;

end
