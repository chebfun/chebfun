% RDIVIDE_LINTEST   Check that correct linearity information is returned for
% various operations with ./

function pass = test_lintest_rdivide()
%% Initialisation
x = chebfun(@(x) x);
u = adchebfun(x) + 2;
v = adchebfun(x) + 3;

%% One variable involved
pass = [];
w = u/2;
pass(length(pass) + 1) = w.linearity == 1;

w = u./(x+2);
pass(length(pass) + 1) = w.linearity == 1;

w = sin(u)./(x+2);
pass(length(pass) + 1) = w.linearity == 0;

w = 2./u;
pass(length(pass) + 1)  = w.linearity == 0;

w = x./u;
pass(length(pass) + 1)  = w.linearity == 0;

w = u./(u+2);
pass(length(pass) + 1)  = w.linearity == 0;

%% Two variables, but still only one in compuations
u = seed(u, 1, 2);
v = seed(v, 2, 2);
w = u/2;
pass(length(pass) + 1)  = all( w.linearity == [1 1]);

w = u./(x+2);
pass(length(pass) + 1)  = all( w.linearity == [1 1]);

w = sin(u)/2;
pass(length(pass) + 1)  = all( w.linearity == [0 1]);

w = 2./u;
pass(length(pass) + 1)  = all( w.linearity == [0 1]);

w = x./u;
pass(length(pass) + 1)  = all( w.linearity == [0 1]);

w = u./(u+2);
pass(length(pass) + 1)  = all( w.linearity == [0 1]);

%% Combination of variables
w = u./v;
pass(length(pass) + 1)  = all( w.linearity == [0 0]);

w = sin(u)./v;
pass(length(pass) + 1)  = all( w.linearity == [0 0]);

w = u./(u.*v);
pass(length(pass) + 1)  = all( w.linearity == [0 0]);

end