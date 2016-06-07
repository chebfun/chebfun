% Test file for cheboppref.m.

function pass = test_cheboppref()

% Test construction from a cheboppref.
p = cheboppref();
pass(1) = isequalNaN(p, cheboppref(p));

% Test construction from a struct.
q = struct();
q.damping = 0;
q.plotting = 'on';
p = cheboppref(q);
pass(2) = (~p.damping) && strcmp(p.plotting, 'on');

% Test overloaded subsasgn.
p = cheboppref();
p.plotting = 'on';
pass(3) = strcmp(p.plotting, 'on');

% Test functions for managing default preferences.
savedPrefs = cheboppref();

cheboppref.setDefaults('factory');
factoryPrefs = cheboppref.getFactoryDefaults();
p = cheboppref();
pass(4) = isequalNaN(p, factoryPrefs);

cheboppref.setDefaults('factory');
p = cheboppref();
p.damping = 0;
p.plotting = 'on';
cheboppref.setDefaults(p);
pass(5) = (~cheboppref().damping) && strcmp(cheboppref().plotting, 'on');

cheboppref.setDefaults('factory');
p = struct();
p.damping = 0;
p.plotting = 'on';
cheboppref.setDefaults(p);
pass(6) = (~cheboppref().damping) && strcmp(cheboppref().plotting, 'on');

cheboppref.setDefaults('factory');
cheboppref.setDefaults('damping', 0, 'plotting', 'on');
pass(7) = (~cheboppref().damping) && strcmp(cheboppref().plotting, 'on');

cheboppref.setDefaults(savedPrefs);

% Test getting defaults:
pass(8) = isnumeric(cheboppref().bvpTol);

% Test specifying discretization via strings
cheboppref.setDefaults('factory');
cheboppref.setDefaults('discretization', 'values');
pref = cheboppref;
pass(9) = strcmp(pref.discretization, 'values');
pref.discretization = 'coeffs';
pass(10) = strcmp(pref.discretization, 'coeffs');


end

function out = isequalNaN(a, b)
    if ( verLessThan('matlab', '7.14') )
        out = isequalwithequalnans(a, b);
    else
        out = isequaln(a, b);
    end
end
