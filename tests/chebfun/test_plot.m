% Tests for chebfun plotting functions.

function pass = test_plot(pref)

if ( nargin < 1 )
    pref = chebpref();
end

% Real scalar functions.
fsr1 = chebfun(@sin, [-1 0 1], pref);
fsr2 = chebfun(@cos, [-1 0.5 1], pref);
fsr3 = chebfun(@exp, [-1 -0.5 1], pref);

% Real array-valued functions.
far1 = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);
far2 = chebfun(@(x) [sin(x) 2*cos(x) 3*exp(x)], [-1 -0.5 1], pref);
far3 = chebfun(@(x) [3*sin(x) 2*cos(x) exp(x)], [-1 0.5 1], pref);

% Real quasimatrices.
fqr1 = cheb2quasi(far1);
fqr2 = cheb2quasi(far2);
fqr3 = cheb2quasi(far3);

% Complex functions.
fsc = chebfun(@(x) exp(1i*x), [-1 0 1], pref);
fac = chebfun(@(x) [exp(1i*x) 2*exp(1i*x) 3*exp(1i*x)], [-1 0.5 1], pref);
fqc = cheb2quasi(fac);

% Singular functions.
fsing = chebfun(@(x) 1./x, [-1 0 1], 'exps', [0 -1 -1 0]);

% Function with a discontinuity.
fdc = chebfun({@sin, @exp}, [-1 0 1], pref);

% Obviously, we can't check if the plots are correct without human
% intervention, so all these tests are meant to do is make sure none of the
% plotting functions crash.

hfig = figure('Visible', 'off');

% Plots of real scalar functions.
pass(1) = doesNotCrash(@() plot(fsr1));
pass(2) = doesNotCrash(@() plot(fsr1, fsr2));
pass(3) = doesNotCrash(@() plot3(fsr1, fsr2, fsr3));

% Plots of real array-valued functions.
pass(4) = doesNotCrash(@() plot(far1));
pass(5) = doesNotCrash(@() plot(far1, far2));
pass(6) = doesNotCrash(@() plot3(far1, far2, far3));

% Plots of real quasimatrices.
pass(7) = doesNotCrash(@() plot(fqr1));
pass(8) = doesNotCrash(@() plot(fqr1, fqr2));
pass(9) = doesNotCrash(@() plot3(fqr1, fqr2, fqr3));

% Plots which mix array-valued functions and quasimatrices.
pass(10) = doesNotCrash(@() plot(far1, fqr2));
pass(11) = doesNotCrash(@() plot(fqr1, far2));

pass(12) = doesNotCrash(@() plot3(far1, far2, fqr3));
pass(13) = doesNotCrash(@() plot3(far1, fqr2, far3));
pass(14) = doesNotCrash(@() plot3(fqr1, far2, far3));
pass(15) = doesNotCrash(@() plot3(fqr1, fqr2, far3));
pass(16) = doesNotCrash(@() plot3(fqr1, far2, fqr3));
pass(17) = doesNotCrash(@() plot3(far1, fqr2, fqr3));

% Plots of quasimatrices against scalar functions.
pass(18) = doesNotCrash(@() plot(fsr1, fqr2));
pass(19) = doesNotCrash(@() plot(fqr1, fsr2));
pass(20) = doesNotCrash(@() plot(fsr1, far2));
pass(21) = doesNotCrash(@() plot(far1, fsr2));

% Plots of complex functions.
pass(22) = doesNotCrash(@() plot(fsc));
pass(23) = doesNotCrash(@() plot(fac));
pass(24) = doesNotCrash(@() plot(fqc));

% Check plot of a singular function.
pass(25) = doesNotCrash(@() plot(fsing));

% Check plot flags and other options.
pass(26) = doesNotCrash(@() plot(fsr1, 'numpts', 100));
pass(27) = doesNotCrash(@() plot(fsr2, 'interval', [-0.5 0.5]));
pass(28) = doesNotCrash(@() plot(fsr2, [-0.5 0.5]));
pass(29) = doesNotCrash(@() plot(fdc, 'jumpline', 'r-'));
pass(30) = doesNotCrash(@() plot(fdc, 'jumpline', 'none'));
pass(31) = doesNotCrash(@() plot3(fdc, fsr1, fsr2, 'jumpline', 'r-'));

% Check plotting discrete data alongside CHEBFUN objects
x = linspace(-1,1,10).';
pass(32) = doesNotCrash(@() plot(far1, 'b', far2, 'r', x, far1(x), 'om', x, far3(x), '-ok'));

% Check SURF, SURFACE, SURFC, and MESH.
pass(33) = doesNotCrash(@() surf(far1));
pass(34) = doesNotCrash(@() surf(fqr1));
% (SURFACE is a wrapper for SURF, so we don't need to be so thorough.)
pass(35) = doesNotCrash(@() surface(fqr1));
pass(36) = doesNotCrash(@() surfc(far1));
pass(37) = doesNotCrash(@() surfc(fqr1));
pass(38) = doesNotCrash(@() mesh(far1));
pass(39) = doesNotCrash(@() mesh(fqr1));

close(hfig);

end

function pass = doesNotCrash(fn)

try
    fn();
    pass = true;
catch ME;
    pass = false;
end

end
