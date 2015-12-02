function pass = test_carrier_US(pref)
% Test Carrier equation.
%
% Toby Driscoll / Asgeir Birkisson, June 2014.
if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e2*pref.bvpTol;
dom = [-1 1];

%%

pref.discretization = @ultraS;
N = chebop(@(x,u) 0.01.*diff(u,2)+2.*(1-x.^2).*u+u.^2-1, dom);
rhs = 0;
N.bc = @(x,u) [u(-1); u(1)];
x = chebfun(@(x) x, dom);
N.init =  2.*(x.^2-1).*(1-2./(1+20.*x.^2));

u = solvebvp(N, rhs, pref);

xx = (-1:.25:1)';
hiquality_ans = [
    0
    -1.487429807540814
    -1.785617248281071
    1.572366197526305
    -1.539652044363185
    1.572366197526230
    -1.785617248281089
    -1.487429807540795
    0
    ];

pass(1) = norm(u(xx) - hiquality_ans) < tol;
pass(2) = norm(u([-1 1])) < tol ;


end


