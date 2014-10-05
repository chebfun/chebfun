% Test file for chebtech1/coeffs2vals.m

function pass = test_coeffs2vals(varargin)

% Set a tolerance (pref.eps doesn't matter)
tol = 100*eps;

%%
% Test that a single coefficient is converted correctly
c = sqrt(2);
v = chebtech1.coeffs2vals(c);
pass(1) = (c == v);

%%
% Simple data (even case)
c = (6:-1:1).';
% Exact values
vTrue = [ -3*sqrt(6)/2-5/sqrt(2)+2*sqrt(3)+7 ; 4 - sqrt(2)/2 ; -3*sqrt(6)/2+5/sqrt(2)-2*sqrt(3)+7 ; 3*sqrt(6)/2-5/sqrt(2)-2*sqrt(3)+7 ; 4 + sqrt(2)/2 ; 3*sqrt(6)/2+5/sqrt(2)+2*sqrt(3)+7];

%%
% Test real branch
v = chebtech1.coeffs2vals(c);
pass(2) = norm(v - vTrue, inf) < tol;
pass(3) = ~any(imag(v));

%%
% Test imaginary branch
v = chebtech1.coeffs2vals(1i*c);
pass(4) = norm(v - 1i*vTrue, inf) < tol;
pass(5) = ~any(real(v));

%%
% Test general branch
v = chebtech1.coeffs2vals((1+1i)*c);
pass(6) = norm(v - (1+1i)*vTrue, inf) < tol;

%%
% Test for array input
v = chebtech1.coeffs2vals([c, -c]);
pass(7) = norm(v(:,1) - vTrue, inf) < tol && ...
          norm(v(:,2) + vTrue, inf) < tol;
      
%%
% Simple data (odd case)
c = (5:-1:1).';
% Exact values
vTrue = [ 11/2+sqrt(5)-2*sqrt((5+sqrt(5))/2)-sqrt((5-sqrt(5))/2) ; 11/2-sqrt(5)-2*sqrt((5-sqrt(5))/2)+sqrt((5+sqrt(5))/2) ; 3 ; 11/2-sqrt(5)+2*sqrt((5-sqrt(5))/2)-sqrt((5+sqrt(5))/2) ; 11/2+sqrt(5)+2*sqrt((5+sqrt(5))/2)+sqrt((5-sqrt(5))/2) ];

%%
% Test real branch
v = chebtech1.coeffs2vals(c);
pass(8) = norm(v - vTrue, inf) < tol;
pass(9) = ~any(imag(v));

%%
% Test imaginary branch
v = chebtech1.coeffs2vals(1i*c);
pass(10) = norm(v - 1i*vTrue, inf) < tol;
pass(11) = ~any(real(v));

%%
% Test general branch
v = chebtech1.coeffs2vals((1+1i)*c);
pass(12) = norm(v - (1+1i)*vTrue, inf) < tol;

%%
% Test for array input
v = chebtech1.coeffs2vals([c, -c]);
pass(13) = norm(v(:,1) - vTrue, inf) < tol && ...
          norm(v(:,2) + vTrue, inf) < tol;
      
end
