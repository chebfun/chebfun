function pass = test_deal(pref)

%% Setup.
if ( nargin == 0 )
    pref = chebfunpref();
end

% Create some CHEBFUNs and OPERATORBLOCKs.
x = chebfun(@(x) x);
f = x;
g = sin(x);
D = operatorBlock.diff;
LD = linop(D);

%% Deal a square CHEBMATRIX.
A = [D, f; g', 2];

a = deal(A);
La = linop(a);

% Check to make sure single output works. Do so by instantiating the converting
% the OPERATORBLOCK to a LINOP, then instantiating it on a grid.
prefs = cheboppref();
prefs.discretization = @chebcolloc2;
pass(1) = norm(matrix(LD, 6, prefs) - matrix(La, 6, prefs)) == 0;

%% Check to make sure multiple outputs work (there is a row-major ordering).
[a,b,c,d] = deal(A);

La = linop(a);
err = norm(matrix(LD, 6, prefs) - matrix(La, 6, prefs)) + norm(f-b) + norm(g'-c) + norm(d-2);
pass(2) = err == 0;

%% Deal a tall CHEBMATRIX.
B = [f; 2*f; 3*f; 4*f];
[a, b, c, d] = deal(B);
abcd = [a;b;c;d];
pass(3) = norm(B-abcd) == 0;


%% Make sure the error works.
pass(4) = 0;
try
    [a,b,c,d,e] = deal(A);
catch ME
    if (ME.identifier == 'CHEBFUN:CHEBMATRIX:deal:tooManyOutputs')
        pass(4) = 1;
    end
end

end
