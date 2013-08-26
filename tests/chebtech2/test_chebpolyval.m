% Test file for chebtech2/chebpolyval.m

function pass = test_chebpolyval(varargin)

% Set a tolerance (pref.eps doesn't matter)
tol = 100*eps;

pass = zeros(1, 7); % Pre-allocate pass matrix.

%%
% Test that a single coefficient is converted correctly
c = sqrt(2);
v = chebtech2.chebpolyval(c);
pass(1) = (c == v);

%%
% Some simple data 
c = (1:5).';
% Exact values
vTrue = [3 ; 4 - sqrt(2) ; 3 ; 4 + sqrt(2) ; 15];

%%
% Test real branch
v = chebtech2.chebpolyval(c);
pass(2) = norm(v - vTrue, inf) < tol;
pass(3) = ~any(imag(v));

%%
% Test imaginary branch
v = chebtech2.chebpolyval(1i*c);
pass(4) = norm(v - 1i*vTrue, inf) < tol;
pass(5) = ~any(real(v));

%%
% Test general branch
v = chebtech2.chebpolyval((1+1i)*c);
pass(6) = norm(v - (1+1i)*vTrue, inf) < tol;

%%
% Test for array input
v = chebtech2.chebpolyval([c, c(end:-1:1)]);
tmp = ones(size(vTrue));
tmp(end-1:-2:1) = -1;
pass(7) = norm(v(:,1) - vTrue, inf) < tol && ...
          norm(v(:,2) - tmp.*vTrue, inf) < tol;
      
end
