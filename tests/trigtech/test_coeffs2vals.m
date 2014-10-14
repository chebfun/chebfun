% Test file for trigtech/coeffs2vals.m

function pass = test_coeffs2vals(varargin)

% Set a tolerance (pref.eps doesn't matter)
tol = 100*eps;

%%
% Test that a single coefficient is converted correctly
c = sqrt(2);
v = trigtech.coeffs2vals(c);
pass(1) = (c == v);

%%
% Simple data (even case)
c = [0 0.5 0 1 0 0.5].';
% Exact values
f_exact = @(x) 1 + cos(2*pi*x);
x_pts = trigtech.trigpts(numel(c));
vTrue = f_exact(x_pts);

%%
% Test real branch
v = trigtech.coeffs2vals(c);
pass(2) = norm(v - vTrue, inf) < tol;
pass(3) = ~any(imag(v));

%%
% Test imaginary branch
v = trigtech.coeffs2vals(1i*c);
pass(4) = norm(v - 1i*vTrue, inf) < tol;
pass(5) = ~any(real(v));

%%
% Test general branch
v = trigtech.coeffs2vals((1+1i)*c);
pass(6) = norm(v - (1+1i)*vTrue, inf) < tol;

%%
% Test for array input
v = trigtech.coeffs2vals([c, -c]);
pass(7) = norm(v(:,1) - vTrue, inf) < tol && ...
          norm(v(:,2) + vTrue, inf) < tol;
      
%%
% Simple data (odd case)
c = [.5 0.5i 0 1 0 -0.5i .5].';
% Exact values
f_exact = @(x) 1 + sin(2*pi*x) + cos(3*pi*x);
x_pts = trigtech.trigpts(numel(c));
vTrue = f_exact(x_pts);

%%
% Test real branch
v = trigtech.coeffs2vals(c);
pass(8) = norm(v - vTrue, inf) < tol;
pass(9) = ~any(imag(v));

%%
% Test imaginary branch
v = trigtech.coeffs2vals(1i*c);
pass(10) = norm(v - 1i*vTrue, inf) < tol;
pass(11) = ~any(real(v));

%%
% Test general branch
v = trigtech.coeffs2vals((1+1i)*c);
pass(12) = norm(v - (1+1i)*vTrue, inf) < tol;

%%
% Test for array input
v = trigtech.coeffs2vals([c, -c]);
pass(13) = norm(v(:,1) - vTrue, inf) < tol && ...
          norm(v(:,2) + vTrue, inf) < tol;
      
end
