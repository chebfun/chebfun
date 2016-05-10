% Test file for SPINPREF/SPINPREF:

function pass = test_spinpref()

% Construction from STRING for AC equation:
pref = spinpref('ac');
Nplot = pref.Nplot;
pass(1) = isequal(Nplot, 1024);

% Construction from inputs:
pref = spinpref('Ylim', [0 1], 'dataToPlot', 'real', 'dealias', 'on');
pass(2) = isequal(pref.Ylim, [0 1]);
pass(3) = strcmpi(pref.dataToPlot, 'real');
pass(4) = strcmpi(pref.dealias, 'on');

pref = spinpref('dt', 1e-1, 'dtmin', 1e-10, 'dtmax', 1, 'errTol', 1e-10);
pass(5) = isequal(pref.dt, 1e-1);
pass(6) = isequal(pref.dtmin, 1e-10);
pass(7) = isequal(pref.dtmax, 1);
pass(8) = isequal(pref.errTol, 1e-10);

pref = spinpref('iterPlot', 10, 'M', 100, 'N', 1024, 'Nmin', 10);
pass(9) = isequal(pref.iterPlot, 10);
pass(10) = isequal(pref.M, 100);
pass(11) = isequal(pref.N, 1024);
pass(12) = isequal(pref.Nmin, 10);

pref = spinpref('Nmax', 10, 'Nplot', 2, 'plot', 'movie', 'scheme', 'lawson4');
pass(13) = isequal(pref.Nmax, 10);
pass(14) = isequal(pref.Nplot, 2);
pass(15) = strcmpi(pref.plot, 'movie');
pass(16) = strcmpi(pref.scheme, 'lawson4');

end
