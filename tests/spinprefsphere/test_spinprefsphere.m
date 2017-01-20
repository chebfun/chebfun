% Test file for SPINPREFSPHERE:

function pass = test_spinprefsphere()

% Construction from inputs:
pref = spinprefsphere('Clim', [0 10], 'dataplot', 'abs', 'dealias', 'off');
pass(1) = isequal(pref.Clim, [0 10]);
pass(2) = strcmpi(pref.dataplot, 'abs');
pass(3) = strcmpi(pref.dealias, 'off');

pref = spinprefsphere('iterplot', 10);
pass(4) = isequal(pref.iterplot, 10);

pref = spinprefsphere('Nplot', 2, 'plot', 'movie');
pass(5) = isequal(pref.Nplot, 2);
pass(6) = strcmpi(pref.plot, 'movie');

pref = spinprefsphere('view', [10 20]);
pass(7) = isequal(pref.view, [10 20]);

end
