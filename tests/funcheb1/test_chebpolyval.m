% Test file for funcheb1/chebpolyval.m

function pass = test_chebpolyval(varargin)

% Set a tolerance (pref.eps doesn't matter)
tol = 100*eps;

%%
% Test that a single coefficient is converted correctly
c = sqrt(2);
v = funcheb1.chebpolyval(c);
pass(1) = (c == v);

%%
% Some simple data 
c = (1:6).';
% Exact values
vTrue = [ -3*sqrt(6)/2-5/sqrt(2)+2*sqrt(3)+7 ; 4 - sqrt(2)/2 ; -3*sqrt(6)/2+5/sqrt(2)-2*sqrt(3)+7 ; 3*sqrt(6)/2-5/sqrt(2)-2*sqrt(3)+7 ; 4 + sqrt(2)/2 ; 3*sqrt(6)/2+5/sqrt(2)+2*sqrt(3)+7];

%%
% Test real branch
v = funcheb1.chebpolyval(c);
pass(2) = norm(v - vTrue, inf) < tol;
pass(3) = ~any(imag(v));

%%
% Test imaginary branch
v = funcheb1.chebpolyval(1i*c);
pass(4) = norm(v - 1i*vTrue, inf) < tol;
pass(5) = ~any(real(v));

%%
% Test general branch
v = funcheb1.chebpolyval((1+1i)*c);
pass(6) = norm(v - (1+1i)*vTrue, inf) < tol;

%%
% Test for array input
v = funcheb1.chebpolyval([c, -c]);
pass(7) = norm(v(:,1) - vTrue, inf) < tol && ...
          norm(v(:,2) + vTrue, inf) < tol;
      
end
