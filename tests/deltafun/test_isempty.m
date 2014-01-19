% Test file for @deltafun/isempty.m

%function pass = test_isempty(pref)

% if (nargin < 1)
%     pref = chebpref();
% end

d = deltafun();
pass(1) = isempty(d);

d = deltafun([], []);
pass(2) = isempty(d);

d = deltafun([], [], []);
pass(3) = isempty(d);

d = deltafun([], [], [], []);
pass(4) = isempty(d);

pass

%end