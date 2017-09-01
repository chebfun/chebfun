function pass = test_plotting()
% Check that the very basic plotting commands do not crash.

% Real-valued functions.
f1 = chebfun3(@(x,y,z) x.*y.*z, [-1 2 -1 2 -1 2]);

% A complex-valued function.
f2 = chebfun3(@(x,y,z) exp(cos(1i*x.*y.*z)));

% Obviously, we can't check if the plots are correct without human
% intervention, so all these tests are meant to do is make sure none of the
% plotting functions crash.

hfig = figure('Visible', 'off');

%% Real-valued function
% Non-GUI plots
pass(1) = doesNotCrash(@() plot(f1));
pass(2) = doesNotCrash(@() slice(f1, 'noslider'));
pass(3) = doesNotCrash(@() slice(f1, 0.5, -0.3, 0.9));
pass(4) = doesNotCrash(@() isosurface(f1, 'noslider'));
pass(5) = doesNotCrash(@() isosurface(f1, [0.5, -0.6]));
pass(6) = doesNotCrash(@() scan(f1));
pass(6) = doesNotCrash(@() scan(f1, 1));
pass(7) = doesNotCrash(@() scan(f1, 1, 'hold'));
pass(8) = doesNotCrash(@() scan(f1, 3));
% GUI-based plots:
pass(9) = doesNotCrash(@() isosurface(f1));
pass(10) = doesNotCrash(@() slice(f1));
pass(11) = doesNotCrash(@() surf(f1));

%% Complex-valued function
% Non-GUI plots
pass(12) = doesNotCrash(@() plot(f2));
pass(13) = doesNotCrash(@() slice(f2, 'noslider'));
pass(14) = doesNotCrash(@() slice(f2, 0.5, -0.3, 0.9));
pass(15) = doesNotCrash(@() isosurface(f2, 'noslider'));
pass(16) = doesNotCrash(@() isosurface(f2, [0.5, -0.6]));
pass(17) = doesNotCrash(@() scan(f2));
pass(18) = doesNotCrash(@() scan(f2, 'hold'));
pass(19) = doesNotCrash(@() scan(f2, 2, 'hold'));
% GUI-based plots:
pass(20) = doesNotCrash(@() isosurface(f2));
pass(21) = doesNotCrash(@() slice(f2));
pass(22) = doesNotCrash(@() surf(f2));

%% Test plotting on a non-standard domain
r = chebfun3(@(r,t,p) r, [0 1 0 2*pi pi/4 pi/2]);
t = chebfun3(@(r,t,p) t, [0 1 0 2*pi pi/4 pi/2]);
p = chebfun3(@(r,t,p) p, [0 1 0 2*pi pi/4 pi/2]);
x = r.*cos(t).*cos(p);
y = r.*sin(t).*cos(p);
z = r.*sin(p);
density = sin(10*t).*cos(10*r)+1;
pass(23) = doesNotCrash(@() plot(x,y,z,density));

close(hfig);

end

function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME
    pass = false;
end

end