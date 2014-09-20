% Test file for chebtech1/vals2coeffs.m

function pass = test_vals2coeffs(varargin)

% Set a tolerance (pref.eps doesn't matter)
tol = 100*eps;

%%
% Test that a single value is converted correctly
v = sqrt(2);
c = chebtech1.vals2coeffs(v);
pass(1) = (v == c);

%%
% Simple data (even case)
v = (1:6).';
% Exact coefficients
cTrue = [ 7/2 ; sqrt(6)/2+5*sqrt(2)/6 ; 0 ; sqrt(2)/6 ; 0 ; sqrt(6)/2-5*sqrt(2)/6 ];
%%
% Test real branch
c = chebtech1.vals2coeffs(v);
pass(2) = norm(c - cTrue, inf) < tol;
pass(3) = ~any(imag(c));

%%
% Test imaginary branch
c = chebtech1.vals2coeffs(1i*v);
pass(4) = norm(c - 1i*cTrue, inf) < tol;
pass(5) = ~any(real(c));

%%
% Test general branch
c = chebtech1.vals2coeffs((1+1i)*v);
pass(6) = norm(c - (1+1i)*cTrue, inf) < tol;

%%
% Test for array input
c = chebtech1.vals2coeffs([v, v(end:-1:1)]);
tmp = ones(size(cTrue));
tmp(2:2:end) = -1;
pass(7) = norm(c(:,1) - cTrue, inf) < tol && ...
          norm(c(:,2) - tmp.*cTrue, inf) < tol;
      
%%
% Simple data (odd case)
v = (1:5).';
% Exact coefficients
cTrue = [3 ; (2/5)*(sqrt((5-sqrt(5))/2)+2*sqrt((5+sqrt(5))/2)); 0 ; (2/5)*(2*sqrt((5-sqrt(5))/2)-sqrt((5+sqrt(5))/2)) ; 0 ;];
%%
% Test real branch
c = chebtech1.vals2coeffs(v);
pass(8) = norm(c - cTrue, inf) < tol;
pass(9) = ~any(imag(c));

%%
% Test imaginary branch
c = chebtech1.vals2coeffs(1i*v);
pass(10) = norm(c - 1i*cTrue, inf) < tol;
pass(11) = ~any(real(c));

%%
% Test general branch
c = chebtech1.vals2coeffs((1+1i)*v);
pass(12) = norm(c - (1+1i)*cTrue, inf) < tol;

%%
% Test for array input
c = chebtech1.vals2coeffs([v, v(end:-1:1)]);
tmp = ones(size(cTrue)); tmp(end-1:-2:1) = -1;
pass(13) = norm(c(:,1) - cTrue, inf) < tol && ...
          norm(c(:,2) - tmp.*cTrue, inf) < tol;
      
end
