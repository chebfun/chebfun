% Test for @chebfun/get.m.

function pass = test_get(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

fs1 = chebfun(@(x) sin(x), [-1 1]);
fs2 = chebfun(@(x) cos(x), [-1 0 1]);

fa1 = chebfun(@(x) [sin(x) cos(x)], [-1 1]);
fa2 = chebfun(@(x) [sin(x) cos(x)], [-1 0 1]);

fq1 = cheb2quasi(fa1);
fq2 = cheb2quasi(fa2);
fq3 = [fs1 fs2];

% Test SIMPLEVEL = 2.
c = get(fs1, 'coeffs');
pass(1) = isnumeric(c) && (size(c, 2) == 1);
c = get(fs2, 'coeffs');
pass(2) = iscell(c) && isequal(size(c), [2 1]);
c = get(fa1, 'coeffs');
pass(3) = isnumeric(c) && (size(c, 2) == 2);
c = get(fa2, 'coeffs');
pass(4) = iscell(c) && isequal(size(c), [2 2]);
c = get(fq1, 'coeffs');
pass(5) = isnumeric(c) && (size(c, 2) == 2);
c = get(fq2, 'coeffs');
pass(6) = iscell(c) && isequal(size(c), [2 2]);
c = get(fq3, 'coeffs');
pass(7) = iscell(c) && isequal(size(c), [1 2]) && ...
    iscell(c{1}) && isequal(size(c{1}), [1 1]) && ...
    iscell(c{2}) && isequal(size(c{2}), [2 1]);

% Test SIMPLEVEL = 1.
c = get(fs1, 'coeffs', 1);
pass(8) = iscell(c) && isequal(size(c), [1 1]);
c = get(fs2, 'coeffs', 1);
pass(9) = iscell(c) && isequal(size(c), [2 1]);
c = get(fa1, 'coeffs', 1);
pass(10) = iscell(c) && isequal(size(c), [1 2]);
c = get(fa2, 'coeffs', 1);
pass(11) = iscell(c) && isequal(size(c), [2 2]);
c = get(fq1, 'coeffs', 1);
pass(12) = iscell(c) && isequal(size(c), [1 2]);
c = get(fq2, 'coeffs', 1);
pass(13) = iscell(c) && isequal(size(c), [2 2]);
c = get(fq3, 'coeffs', 1);
pass(14) = iscell(c) && isequal(size(c), [1 2]) && ...
    iscell(c{1}) && isequal(size(c{1}), [1 1]) && ...
    iscell(c{2}) && isequal(size(c{2}), [2 1]);

% Test SIMPLEVEL = 0.
c = get(fs1, 'coeffs', 0);
pass(15) = iscell(c) && isequal(size(c), [1 1]) && ...
    iscell(c{1}) && isequal(size(c{1}), [1 1]);
c = get(fs2, 'coeffs', 0);
pass(16) = iscell(c) && isequal(size(c), [1 1]) && ...
    iscell(c{1}) && isequal(size(c{1}), [2 1]);
c = get(fa1, 'coeffs', 0);
pass(17) = iscell(c) && isequal(size(c), [1 2]) && ...
    iscell(c{1}) && isequal(size(c{1}), [1 1]) && ...
    iscell(c{2}) && isequal(size(c{2}), [1 1]);
c = get(fa2, 'coeffs', 0);
pass(18) = iscell(c) && isequal(size(c), [1 2]) && ...
    iscell(c{1}) && isequal(size(c{1}), [2 1]) && ...
    iscell(c{2}) && isequal(size(c{2}), [2 1]);
c = get(fq1, 'coeffs', 0);
pass(19) = iscell(c) && isequal(size(c), [1 2]) && ...
    iscell(c{1}) && isequal(size(c{1}), [1 1]) && ...
    iscell(c{2}) && isequal(size(c{2}), [1 1]);
c = get(fq2, 'coeffs', 0);
pass(20) = iscell(c) && isequal(size(c), [1 2]) && ...
    iscell(c{1}) && isequal(size(c{1}), [2 1]) && ...
    iscell(c{2}) && isequal(size(c{2}), [2 1]);
c = get(fq3, 'coeffs', 0);
pass(21) = iscell(c) && isequal(size(c), [1 2]) && ...
    iscell(c{1}) && isequal(size(c{1}), [1 1]) && ...
    iscell(c{2}) && isequal(size(c{2}), [2 1]);

% Check behavior for row CHEBFUNs.  (SIMPLEVEL = 2 only.)
c = get(fs1.', 'coeffs');
pass(22) = isnumeric(c) && (size(c, 1) == 1);
c = get(fs2.', 'coeffs');
pass(23) = iscell(c) && isequal(size(c), [1 2]);
c = get(fa1.', 'coeffs');
pass(24) = isnumeric(c) && (size(c, 1) == 2);
c = get(fa2.', 'coeffs');
pass(25) = iscell(c) && isequal(size(c), [2 2]);
c = get(fq1.', 'coeffs');
pass(26) = isnumeric(c) && (size(c, 1) == 2);
c = get(fq2.', 'coeffs');
pass(27) = iscell(c) && isequal(size(c), [2 2]);
c = get(fq3.', 'coeffs');
pass(28) = iscell(c) && isequal(size(c), [2 1]) && ...
    iscell(c{1}) && isequal(size(c{1}), [1 1]) && ...
    iscell(c{2}) && isequal(size(c{2}), [2 1]);

% Check getting of delta functions.
x = chebfun(@(x) x);
fd1 = fs1 + dirac(x);
fd2 = fs2 + dirac(x - 0.5) + dirac(x + 0.5);
fd3 = [fd1 fd2];

c = get(fd1, 'deltas');
pass(29) = ( norm(c - [0 ; 1], inf) < 5*eps );
c = get(fd2, 'deltas');
pass(30) = ( norm(c - [-0.5 0.5 ; 1 1], inf) < 5*eps );
c = get(fd3, 'deltas');
pass(31) = iscell(c) && isequal(size(c), [1 2]) && ...
    ( norm(c{1} - [0 ; 1], inf) < 5*eps ) && ...
    ( norm(c{2} - [-0.5 0.5 ; 1 1], inf) < 5*eps );

% Check getting of exponents.
fse1 = chebfun(@(x) 1./(x + 2), [-2 2], 'exps', [-1 0]);
fse2 = chebfun(@(x) gamma(x), [-2 -1 0 2], 'exps', [-1 -1 -1 -1 -1 0]);
fqe1 = [fse2 fse2];
fqe2 = [fse1 fse2];

exps = get(fse1, 'exponents');
pass(32) = isequal(exps, [-1 0]);
exps = get(fse1, 'exponents', 1);
pass(33) = iscell(exps) && isequal(size(exps), [1 1]);
exps = get(fse1, 'exponents', 0);
pass(34) = iscell(exps) && isequal(size(exps), [1 1]) && ...
    iscell(exps{1}) && isequal(size(exps{1}), [1 1]);
exps = get(fse2, 'exponents');
pass(35) = isequal(exps, [-1 -1 ; -1 -1 ; -1 0]);
exps = get(fse2, 'exponents', 1);
pass(36) = iscell(exps) && isequal(size(exps), [3 1]);
exps = get(fqe1, 'exponents');
pass(37) = iscell(exps) && isequal(size(exps), [3 2]);
exps = get(fqe2, 'exponents');
pass(38) = iscell(exps) && isequal(size(exps), [1 2]) && ...
    iscell(exps{1}) && isequal(size(exps{1}), [1 1]) && ...
    iscell(exps{2}) && isequal(size(exps{2}), [3 1]);

end
