% Test file for chebtech2/coeffs2vals.m

function pass = test_coeffs2vals(varargin)

% Set a tolerance (pref.chebfuneps doesn't matter)
tol = 100*eps;

%%
% Test that a single coefficient is converted correctly
c = sqrt(2);
v = chebtech2.coeffs2vals(c);
pass(1) = (c == v);

%%
% Some simple data 
c = (1:5).';
% Exact values
vTrue = [ 3; -4 + sqrt(2); 3; -4-sqrt(2); 15];
%%
% Test real branch
v = chebtech2.coeffs2vals(c);
pass(2) = norm(v - vTrue, inf) < tol;
pass(3) = ~any(imag(v));

%%
% Test imaginary branch
v = chebtech2.coeffs2vals(1i*c);
pass(4) = norm(v - 1i*vTrue, inf) < tol;
pass(5) = ~any(real(v));

%%
% Test general branch
v = chebtech2.coeffs2vals((1+1i)*c);
pass(6) = norm(v - (1+1i)*vTrue, inf) < tol;

%%
% Test for array input
v = chebtech2.coeffs2vals([c, c(end:-1:1)]);
tmp = ones(size(vTrue));
tmp(end-1:-2:1) = -1;
pass(7) = norm(v(:,1) - vTrue, inf) < tol && ...
          norm(v(:,2) - tmp.*vTrue, inf) < tol;
      
%%
% Test for symmetry preservation
c = kron(ones(10,1),eye(2));
v = chebtech1.coeffs2vals(c);
pass(8) = norm(v(:,1) - flipud(v(:,1)), inf) == 0 && ...
          norm(v(:,2) + flipud(v(:,2)), inf) == 0;
      
end
