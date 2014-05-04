% Test file for cheboppref.m.

function pass = test_cheboppref()

% Test construction from a cheboppref.
p = cheboppref();
pass(1) = isequalNaN(p, cheboppref(p));

% Test construction from a struct.
q = struct();
q.damped = 0;
q.plotting = 'on';
p = cheboppref(q);
pass(2) = (~p.damped) && strcmp(p.plotting, 'on');

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
p.damped = 0;
p.plotting = 'on';
cheboppref.setDefaults(p);
pass(5) = (~cheboppref().damped) && strcmp(cheboppref().plotting, 'on');

cheboppref.setDefaults('factory');
p = struct();
p.damped = 0;
p.plotting = 'on';
cheboppref.setDefaults(p);
pass(6) = (~cheboppref().damped) && strcmp(cheboppref().plotting, 'on');

cheboppref.setDefaults('factory');
cheboppref.setDefaults('damped', 0, 'plotting', 'on');
pass(7) = (~cheboppref().damped) && strcmp(cheboppref().plotting, 'on');

cheboppref.setDefaults(savedPrefs);

end

function out = isequalNaN(a, b)
    if ( verLessThan('matlab', '7.14') )
        out = isequalwithequalnans(a, b);
    else
        out = isequaln(a, b);
    end
end
