% Test file for SPINPREF2:

function pass = test_spinpref2()

% Construction from inputs:
pref = spinpref2('Clim', [0 10], 'dataplot', 'abs', 'dealias', 'off');
pass(1) = isequal(pref.Clim, [0 10]);
pass(2) = strcmpi(pref.dataplot, 'abs');
pass(3) = strcmpi(pref.dealias, 'off');

pref = spinpref2('iterplot', 10, 'M', 100);
pass(4) = isequal(pref.iterplot, 10);
pass(5) = isequal(pref.M, 100);

pref = spinpref2('Nplot', 2, 'plot', 'movie', 'scheme', 'lawson4');
pass(6) = isequal(pref.Nplot, 2);
pass(7) = strcmpi(pref.plot, 'movie');
pass(8) = strcmpi(pref.scheme, 'lawson4');

pref = spinpref2('view', [10 20]);
pass(9) = isequal(pref.view, [10 20]);

end
