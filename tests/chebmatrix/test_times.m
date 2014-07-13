function pass = test_times(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% This test comers from #806.
x = chebfun(@(x) x, pref);
xMat = [x, 2*x ; 3*x, 4*x];
A = [1 2 ; 3 4]*xMat;
B = [7*x 10*x ; 15*x 22*x];
pass(1) = norm(A - B) < 1e-10;

%% Test from #1013.
f = chebmatrix([sin(x), cos(x)]);
g = [2,3].*f; 
pass(2) = numColumns(g{1}) == 1;
pass(3) = norm(g{1} - 2*sin(x)) + norm(g{2} - 3*cos(x)) == 0;

% Make sure addition of the CHEBMATRIX objects is OK.
h = f + g;
pass(4) = norm(h{1} - 3*sin(x)) + norm(h{2} - 4*cos(x)) == 0;
end