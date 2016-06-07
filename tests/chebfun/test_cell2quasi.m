function pass = test_cell2quasi(pref)
%TEST_CELL2QUASI
%
% Test that a quasimatrix created from a chebfun and a scalar has the correct
% domain.
if ( nargin == 0 )
    pref = chebfunpref();
end

% See #1344.
x35 = chebfun(@(x) x, [3, 5], pref);
Q = quasimatrix({x35, 2});
pass(1) = all( Q(:,2).domain == [3, 5] );


end
