% Test file for chebfunpref.m.

function pass = test_chebfunpref()

% Test construction from a chebfunpref.
p = chebfunpref();
pass(1) = isequalNaN(p, chebfunpref(p));

% Test construction from a struct.
q = struct();
q.enableBreakpointDetection = true;
q.testPref = 'test';
p = chebfunpref(q);
pass(2) = p.enableBreakpointDetection && strcmp(p.techPrefs.testPref, 'test');

% Test construction from a struct in which incomplete techPrefs substructures
% will need to be merged.
q = struct();
q.testPref1 = 'test1';
q.techPrefs.testPref2 = 'test2';
p = chebfunpref(q);
pass(3) = strcmp(p.techPrefs.testPref1, 'test1') && ...
    strcmp(p.techPrefs.testPref2, 'test2');

% Test overloaded subsref.
q = struct();
q.techPrefs.testPref = 'test';
q.techPrefs.subPrefs.testSubPref = 'subTest';
p = chebfunpref(q);
pass(4) = strcmp(p.techPrefs.testPref, 'test') && strcmp(p.testPref, 'test');
pass(5) = strcmp(p.techPrefs.subPrefs.testSubPref, 'subTest') && ...
    strcmp(p.subPrefs.testSubPref, 'subTest');

% Test overloaded subsasgn.
p = chebfunpref();
p.maxLength = 1337;
pass(6) = p.maxLength == 1337;

p.techPrefs.testPref1 = 'test1';
pass(7) = strcmp(p.techPrefs.testPref1, 'test1');

p.testPref2 = 'test2';
pass(8) = strcmp(p.techPrefs.testPref2, 'test2');

% Test behavior of mergePrefs().
p = struct();
p.testPref = 'test';
q = struct();
pass(9) = isequalNaN(chebfunpref.mergePrefs(p, q), p);

q.testPref = 'testq';
pass(10) = strcmp(chebfunpref.mergePrefs(p, q).testPref, 'testq');

p.testPref2 = 'test2';
q.testPref2q = 'test2q';
map.testPref2q = 'testPref2';
pass(11) = strcmp(chebfunpref.mergePrefs(p, q, map).testPref2, 'test2q');

% Test behavior of mergePrefs() for chebfunpref inputs.
p = chebfunpref();
p.techPrefs.testPref = 'test';
q = struct();
q.testPref = 'testq';
pass(12) = isequalNaN(chebfunpref.mergePrefs(p, q), ...
    chebfunpref.mergePrefs(p.techPrefs, q));
pass(13) = isequalNaN(chebfunpref.mergePrefs(q, p), ...
    chebfunpref.mergePrefs(q, p.techPrefs));

q = chebfunpref();
q.techPrefs.testPref = 'testq';
pass(14) = isequalNaN(chebfunpref.mergePrefs(p, q), ...
    chebfunpref.mergePrefs(p.techPrefs, q.techPrefs));

% Test functions for managing default preferences.
savedPrefs = chebfunpref();

chebfunpref.setDefaults('factory');
factoryPrefs = chebfunpref.getFactoryDefaults();
p = chebfunpref();
pass(15) = isequalNaN(p, factoryPrefs);

chebfunpref.setDefaults('factory');
p = chebfunpref();
p.domain = [-2 7];
p.testPref = 'testq';
chebfunpref.setDefaults(p);
pass(16) = strcmp(chebfunpref().testPref, 'testq') && ...
    isequal(chebfunpref().domain, [-2 7]);

chebfunpref.setDefaults('factory');
p = struct();
p.domain = [-2 7];
p.testPref = 'testq';
chebfunpref.setDefaults(p);
pass(17) = strcmp(chebfunpref().testPref, 'testq') && ...
    isequal(chebfunpref().domain, [-2 7]);

chebfunpref.setDefaults('factory');
chebfunpref.setDefaults('domain', [-2 7], 'testPref', 'testq');
pass(18) = strcmp(chebfunpref().testPref, 'testq') && ...
    isequal(chebfunpref().domain, [-2 7]);

chebfunpref.setDefaults(savedPrefs);

end

function out = isequalNaN(a, b)
    if ( verLessThan('matlab', '7.14') )
        out = isequalwithequalnans(a, b);
    else
        out = isequaln(a, b);
    end
end
