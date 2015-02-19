function pass = test_cell2quasi(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% See #1344.
x35 = chebfun(@(x) x, [3, 5], pref);
Q = quasimatrix({x35, 2});
pass(1) = all( Q(:,2).domain == [3, 5] );

% x = chebfun(@(x) x, [2, 4]);
% x18 = [x ; 1.8];
% plot(x18)

end