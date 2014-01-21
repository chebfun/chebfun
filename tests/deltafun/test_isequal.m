% Test file for @deltafun/isequal.m

function pass = test_isequal(pref)

% if (nargin < 1)
%     pref = chebpref();
% end
%%
d1 = deltafun();
d2 = deltafun();

pass(1) = isequal(d1, []) && isequal( [], d1) && isequal( d1, d2); 

f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = rand(1,5);

d1 = deltafun(f, mag, loc);
d2 = deltafun(f, mag, loc);

pass(2) = isequal(d1, d2);

pass(3) = ~isequal( d1, .992312341234*d2 );

d2.impulses = d2.impulses(1:end-1, :);
pass(4) = ~isequal( d1, d2);

d1.impulses = d1.impulses(:, 1:end-1);
d1.location = d1.location(1:end-1);

pass(5) = ~isequal(d1, d2);

d1 = deltafun( [], 1, 0);
d2 = deltafun( [], [1; 0], 0 );
pass(6) = isequal(d1, d2);

d1 = deltafun(f, [], [] );
d2 = deltafun([], [], []);

pass(7) = ~isequal(d1, d2);

end