function pass = test_partitionCombine( ) 
% Test the partition and combine methods.

tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;

% Check that partitioning an empty spherefun gives two empty spherefuns
f = spherefun;
[feven,fodd] = partition(f);
pass(1) = isempty(feven) & isempty(fodd);

feven = spherefun;
fodd = spherefun;
f = combine(feven,fodd);
pass(2) = isempty(f);

% Check that an even spherefun is partitioned correctly.
f = spherefun(@(x,y,z) sin(pi*x.*y));  % Strictly even/pi-periodic
[feven,fodd] = partition(f);
pass(3) = isequal(feven,f);
pass(4) = isempty(fodd);

% Check that odd spherefun is partitioned correctly.
f = spherefun(@(x,y,z) sin(pi*x.*z));  % Strictly odd/anti-periodic
[feven,fodd] = partition(f);
pass(5) = isequal(fodd,f);
pass(6) = isempty(feven);

% Check that the even and odd terms are partitioned correctly
fe = @(x,y,z) sin(pi*x.*y);  % Strictly even/pi-periodic
fo = @(x,y,z) sin(pi*x.*z);  % Strictly odd/anti-periodic
f = spherefun(@(x,y,z) fe(x,y,z) + fo(x,y,z));
[feven,fodd] = partition(f);
pass(7) = norm(spherefun(fe)-feven) < tol;
pass(8) = norm(spherefun(fo)-fodd) < tol;

% Check that combine puts an even spherefun and empty spherfun back
% together.
feven = spherefun(@(x,y,z) sin(pi*x.*y));  % Strictly even/pi-periodic
fodd = spherefun;
f = combine(feven,fodd);
pass(9) = norm(feven-f) < tol;

% Check that combine puts an even spherefun and empty spherfun back
% together.
feven = spherefun;
fodd = spherefun(@(x,y,z) sin(pi*x.*z));  % Strictly odd/anti-periodic
f = combine(feven,fodd);
pass(10) = norm(fodd-f) < tol;

% Check that combine puts non-empty odd and even spherefuns back together.
fe = @(x,y,z) sin(pi*x.*y);  % Strictly even/pi-periodic
fo = @(x,y,z) sin(pi*x.*z);  % Strictly odd/anti-periodic
fcombine = spherefun(@(x,y,z) fe(x,y,z) + fo(x,y,z));
f = combine(spherefun(fe),spherefun(fo));
pass(11) = norm(fcombine-f) < tol;

% Check that combine cannot put together spherefuns that are not strictly
% odd or even.
try
    f = spherefun(@(x,y,z) sin(pi*x.*y) + sin(pi*x.*z));   
    g = combine(f,f);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:SPHEREFUN:combine:parity');
end

end