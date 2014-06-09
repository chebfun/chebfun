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

% Test construction with a non-default base set of preferences and a
% chebfunpref with differences that need to be merged into that set.
p = chebfunpref();
p.domain = [-2 7];
p.eps = 1.0e-6;
p.maxLength = 1337;
p.bogusPref = true;

q = chebfunpref();
q.enableBreakpointDetection = true;
q.eps = 1.0e-12;

r = chebfunpref(p, q);
pass(9) = isequal(r.domain, q.domain) && ...
    isequal(r.enableBreakpointDetection, q.enableBreakpointDetection) && ...
    isequal(r.eps, q.eps) && ...
    isequal(r.maxLength, q.maxLength) && ...
    isequal(r.bogusPref, p.bogusPref);

% Test construction with a non-default base set of preferences and a struct
% with differences that need to be merged into that set.
q = struct();
q.enableBreakpointDetection = true;
q.eps = 1.0e-12;

r = chebfunpref(p, q);
pass(10) = isequal(r.domain, p.domain) && ...
    isequal(r.enableBreakpointDetection, q.enableBreakpointDetection) && ...
    isequal(r.eps, q.eps) && ...
    isequal(r.maxLength, p.maxLength) && ...
    isequal(r.bogusPref, p.bogusPref);

% Test behavior of mergePrefs().
p = struct();
p.testPref = 'test';
q = struct();
pass(11) = isequalNaN(chebfunpref.mergePrefs(p, q), p);

q.testPref = 'testq';
pass(12) = strcmp(chebfunpref.mergePrefs(p, q).testPref, 'testq');

p.testPref2 = 'test2';
q.testPref2q = 'test2q';
map.testPref2q = 'testPref2';
pass(13) = strcmp(chebfunpref.mergePrefs(p, q, map).testPref2, 'test2q');

% Test behavior of mergePrefs() for chebfunpref inputs.
p = chebfunpref();
p.techPrefs.testPref = 'test';
q = struct();
q.testPref = 'testq';
pass(14) = isequalNaN(chebfunpref.mergePrefs(p, q), ...
    chebfunpref.mergePrefs(p.techPrefs, q));
pass(15) = isequalNaN(chebfunpref.mergePrefs(q, p), ...
    chebfunpref.mergePrefs(q, p.techPrefs));

q = chebfunpref();
q.techPrefs.testPref = 'testq';
pass(16) = isequalNaN(chebfunpref.mergePrefs(p, q), ...
    chebfunpref.mergePrefs(p.techPrefs, q.techPrefs));

% Test functions for managing default preferences.
savedPrefs = chebfunpref();

try 
    chebfunpref.setDefaults('factory');
    factoryPrefs = chebfunpref.getFactoryDefaults();
    p = chebfunpref();
    pass(17) = isequalNaN(p, factoryPrefs);

    chebfunpref.setDefaults('factory');
    p = chebfunpref();
    p.domain = [-2 7];
    p.testPref = 'testq';
    chebfunpref.setDefaults(p);
    pass(18) = strcmp(chebfunpref().testPref, 'testq') && ...
        isequal(chebfunpref().domain, [-2 7]);

    chebfunpref.setDefaults('factory');
    p = struct();
    p.domain = [-2 7];
    p.testPref = 'testq';
    chebfunpref.setDefaults(p);
    pass(19) = strcmp(chebfunpref().testPref, 'testq') && ...
        isequal(chebfunpref().domain, [-2 7]);

    chebfunpref.setDefaults('factory');
    chebfunpref.setDefaults('domain', [-2 7], 'testPref', 'testq');
    pass(20) = strcmp(chebfunpref().testPref, 'testq') && ...
        isequal(chebfunpref().domain, [-2 7]);

    % Test getting defaults:
    pass(21) = isnumeric(chebfunpref().eps);
    pass(22) = ischar(chebfunpref().singPrefs.defaultSingType);
    pass(23) = ischar(chebfunpref().refinementFunction);
    
catch ME
    
    % Reset the preferences:
    chebfunpref.setDefaults(savedPrefs);
    
    % Rethrow the error:
    rethrow(ME)
    
end

% Reset the preferences:
chebfunpref.setDefaults(savedPrefs);

end

function out = isequalNaN(a, b)
    if ( verLessThan('matlab', '7.14') )
        out = isequalwithequalnans(a, b);
    else
        out = isequaln(a, b);
    end
end
