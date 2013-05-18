% Test file for bndfun constructor.

function pass = test_createmap(pref)

% Get preferences:
if ( nargin < 1 )
    pref = bndfun.pref;
end
% Set the tolerance:
pref = chebtech.pref(pref);

pass = zeros(1, 2); % Pre-allocate pass matrix.

dom = [-2 7];
map = bndfun.createMap(dom);

pass(1) = ( all( map.for([-1 1]) == dom ) && all( map.inv(dom) == [-1 1] ) )

pass(2) = map.der(1) == diff(dom)/2