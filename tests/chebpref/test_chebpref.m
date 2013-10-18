% Test file for chebpref.m.

function pass = test_chebpref()

% Test construction from a chebpref.
p = chebpref();
pass(1) = isequal(p, chebpref(p));

% Test construction from a struct.
q = struct();
q.enableBreakpointDetection = true;
q.testPref = 'test';
p = chebpref(q);
pass(2) = p.enableBreakpointDetection && strcmp(p.techPrefs.testPref, 'test');

% Test construction from a struct in which incomplete techPrefs substructures
% will need to be merged.
q = struct();
q.testPref1 = 'test1';
q.techPrefs.testPref2 = 'test2';
p = chebpref(q);
pass(3) = isequal(p.techPrefs, ...
    struct('testPref1', 'test1', 'testPref2', 'test2'));

% Test overloaded subsref.
q = struct();
q.techPrefs.testPref = 'test';
q.techPrefs.subPrefs.testSubPref = 'subTest';
p = chebpref(q);
pass(4) = strcmp(p.techPrefs.testPref, 'test') && strcmp(p.testPref, 'test');
pass(5) = strcmp(p.techPrefs.subPrefs.testSubPref, 'subTest') && ...
    strcmp(p.subPrefs.testSubPref, 'subTest');

% Test overloaded subsasgn.
p = chebpref();
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
pass(9) = isequal(chebpref.mergePrefs(p, q), p);

q.testPref = 'testq';
pass(10) = strcmp(chebpref.mergePrefs(p, q).testPref, 'testq');

p.testPref2 = 'test2';
q.testPref2q = 'test2q';
map.testPref2q = 'testPref2';
pass(11) = strcmp(chebpref.mergePrefs(p, q, map).testPref2, 'test2q');

end
