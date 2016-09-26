function pass = test_partitionCombine( ) 
% Test the partition and combine methods.

tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;

% Check that partitioning an empty diskfun gives two empty diskfuns
f = diskfun;
[feven,fodd] = partition(f);
pass(1) = isempty(feven) & isempty(fodd);

feven = diskfun;
fodd = diskfun;
f = combine(feven,fodd);
pass(2) = isempty(f);

% Check that an even diskfun is partitioned correctly.
f = diskfun(@(x,y) sin(pi*x.*y));  % Strictly even/pi-periodic
[feven,fodd] = partition(f);
pass(3) = isequal(feven,f);
pass(4) = isempty(fodd);

% Check that odd diskfun is partitioned correctly.
f = diskfun(@(x,y) sin(pi*x));  % Strictly odd/anti-periodic
[feven,fodd] = partition(f);
pass(5) = isequal(fodd,f);
pass(6) = isempty(feven);

% Check that the even and odd terms are partitioned correctly
fe = @(x,y) sin(pi*x.*y);  % Strictly even/pi-periodic
fo = @(x,y) sin(pi*x);  % Strictly odd/anti-periodic
f = diskfun(@(x,y) fe(x,y) + fo(x,y));
[feven,fodd] = partition(f);
pass(7) = norm(diskfun(fe)-feven) < tol;
pass(8) = norm(diskfun(fo)-fodd) < tol;

% Check that combine puts an even diskfun and empty spherfun back
% together.
feven = diskfun(@(x,y) sin(pi*x.*y));  % Strictly even/pi-periodic
fodd = diskfun;
f = combine(feven,fodd);
pass(9) = norm(feven-f) < tol;

% Check that combine puts an even diskfun and empty spherfun back
% together.
feven = diskfun;
fodd = diskfun(@(x,y) sin(pi*x));  % Strictly odd/anti-periodic
f = combine(feven,fodd);
pass(10) = norm(fodd-f) < tol;

% Check that combine puts non-empty odd and even diskfuns back together.
fe = @(x,y) sin(pi*x.*y);  % Strictly even/pi-periodic
fo = @(x,y) sin(pi*x);  % Strictly odd/anti-periodic
fcombine = diskfun(@(x,y) fe(x,y) + fo(x,y));
f = combine(diskfun(fe),diskfun(fo));
pass(11) = norm(fcombine-f) < tol;

% Check that combine cannot put together diskfuns that are not strictly
% odd or even.
try
    f = diskfun(@(x,y) sin(pi*x.*y) + sin(pi*x));   
    g = combine(f,f);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:DISKFUN:combine:parity');
end

end