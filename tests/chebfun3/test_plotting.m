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
pass(7) = doesNotCrash(@() scan(f1, 'hold'));
pass(8) = doesNotCrash(@() scan(f1, 1));
pass(9) = doesNotCrash(@() scan(f1, 1, 'hold'));
pass(10) = doesNotCrash(@() scan(f1, 2));
pass(11) = doesNotCrash(@() scan(f1, 3));
% GUI-based plots:
pass(12) = doesNotCrash(@() isosurface(f1));
pass(13) = doesNotCrash(@() slice(f1));
pass(14) = doesNotCrash(@() surf(f1));

%% Complex-valued function
% Non-GUI plots
pass(15) = doesNotCrash(@() plot(f2));
pass(16) = doesNotCrash(@() slice(f2, 'noslider'));
pass(17) = doesNotCrash(@() slice(f2, 0.5, -0.3, 0.9));
pass(18) = doesNotCrash(@() isosurface(f2, 'noslider'));
pass(19) = doesNotCrash(@() isosurface(f2, [0.5, -0.6]));
pass(20) = doesNotCrash(@() scan(f2));
pass(21) = doesNotCrash(@() scan(f2, 'hold'));
pass(22) = doesNotCrash(@() scan(f2, 1));
pass(23) = doesNotCrash(@() scan(f2, 1, 'hold'));
pass(24) = doesNotCrash(@() scan(f2, 2, 'hold'));
pass(25) = doesNotCrash(@() scan(f2, 3, 'hold'));
% GUI-based plots:
pass(26) = doesNotCrash(@() isosurface(f2));
pass(27) = doesNotCrash(@() slice(f2));
pass(28) = doesNotCrash(@() surf(f2));

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