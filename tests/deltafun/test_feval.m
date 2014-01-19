% Test file for @deltafun/feval.m

%function pass = test_feval(pref)

% if ( nargin < 1 )
%     pref = ?
% end
d = deltafun();
pass(1) = isempty(d);

x 

d = deltafun([], []);
pass(2) = isempty(d);

d = deltafun([], [], []);
pass(3) = isempty(d);

d = deltafun([], [], [], []);
pass(4) = isempty(d);

pass

pass
%end