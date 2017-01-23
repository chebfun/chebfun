% Test file for SPINPREF:

function pass = test_spinpref()

% Construction from inputs:
pref = spinpref('Ylim', [0 1], 'dataplot', 'real', 'dealias', 'on');
pass(1) = isequal(pref.Ylim, [0 1]);
pass(2) = strcmpi(pref.dataplot, 'real');
pass(3) = strcmpi(pref.dealias, 'on');

pref = spinpref('iterplot', 10, 'M', 100);
pass(4) = isequal(pref.iterplot, 10);
pass(5) = isequal(pref.M, 100);

pref = spinpref('Nplot', 2, 'plot', 'movie', 'scheme', 'lawson4');
pass(6) = isequal(pref.Nplot, 2);
pass(7) = strcmpi(pref.plot, 'movie');
pass(8) = strcmpi(pref.scheme, 'lawson4');

end
