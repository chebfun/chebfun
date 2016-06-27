function pass = test_plotting()
% Check that the very basic plotting commands do not crash.

% Real-valued functions.
f1 = chebfun3(@(x,y,z) exp(cos(10*x.*y.*z)));
f2 = chebfun3(@(x,y,z) x.*y.*z, [-1 2 -1 2 -1 2]);

% A complex-valued function.
f3 = chebfun3(@(x,y,z) exp(cos(1i*x.*y.*z)));

% Obviously, we can't check if the plots are correct without human
% intervention, so all these tests are meant to do is make sure none of the
% plotting functions crash.

hfig = figure('Visible', 'off');

% % Non-GUI plots
pass(1) = doesNotCrash(@() plot(f1));
pass(2) = doesNotCrash(@() slice(f1, 'noslider'));
pass(3) = doesNotCrash(@() slice(f1, 0.5, -0.3, 0.9));
pass(4) = doesNotCrash(@() isosurface(f1, 'noslider'));
pass(5) = doesNotCrash(@() isosurface(f1, [0.5, -0.6]));
pass(6) = doesNotCrash(@() scan(f1));
pass(7) = doesNotCrash(@() scan(f1, 'hold'));
pass(8) = doesNotCrash(@() scan(f1, 1));
pass(9) = doesNotCrash(@() scan(f1, 1, 'hold'));
pass(10) = doesNotCrash(@() scan(f1, 2));
pass(11) = doesNotCrash(@() scan(f1, 2, 'hold'));
pass(12) = doesNotCrash(@() scan(f1, 3));
pass(13) = doesNotCrash(@() scan(f1, 3, 'hold'));
% GUI-based plots:
pass(14) = doesNotCrash(@() isosurface(f1));
pass(15) = doesNotCrash(@() slice(f1));
pass(16) = doesNotCrash(@() surf(f1));


% % Non-GUI plots
pass(17) = doesNotCrash(@() plot(f2));
pass(18) = doesNotCrash(@() slice(f2, 'noslider'));
pass(19) = doesNotCrash(@() slice(f2, 0.5, -0.3, 0.9));
pass(20) = doesNotCrash(@() isosurface(f2, 'noslider'));
pass(21) = doesNotCrash(@() isosurface(f2, [0.5, -0.6]));
pass(22) = doesNotCrash(@() scan(f2));
pass(23) = doesNotCrash(@() scan(f2, 'hold'));
pass(24) = doesNotCrash(@() scan(f2, 1));
pass(25) = doesNotCrash(@() scan(f2, 1, 'hold'));
pass(26) = doesNotCrash(@() scan(f2, 2));
pass(27) = doesNotCrash(@() scan(f2, 2, 'hold'));
pass(28) = doesNotCrash(@() scan(f2, 3));
pass(29) = doesNotCrash(@() scan(f2, 3, 'hold'));
% GUI-based plots:
pass(30) = doesNotCrash(@() isosurface(f2));
pass(31) = doesNotCrash(@() slice(f2));
pass(32) = doesNotCrash(@() surf(f2));


% % Non-GUI plots
pass(33) = doesNotCrash(@() plot(f3));
pass(34) = doesNotCrash(@() slice(f3, 'noslider'));
pass(35) = doesNotCrash(@() slice(f3, 0.5, -0.3, 0.9));
pass(36) = doesNotCrash(@() isosurface(f3, 'noslider'));
pass(37) = doesNotCrash(@() isosurface(f3, [0.5, -0.6]));
pass(38) = doesNotCrash(@() scan(f3));
pass(39) = doesNotCrash(@() scan(f3, 'hold'));
pass(40) = doesNotCrash(@() scan(f3, 1));
pass(41) = doesNotCrash(@() scan(f3, 1, 'hold'));
pass(42) = doesNotCrash(@() scan(f3, 2));
pass(43) = doesNotCrash(@() scan(f3, 2, 'hold'));
pass(44) = doesNotCrash(@() scan(f3, 3));
pass(45) = doesNotCrash(@() scan(f3, 3, 'hold'));
% GUI-based plots:
pass(46) = doesNotCrash(@() isosurface(f3));
pass(47) = doesNotCrash(@() slice(f3));
pass(48) = doesNotCrash(@() surf(f3));

close(hfig);

end


function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end

end