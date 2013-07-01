% Test file for bndfun/createMap.m.

function pass = test_createMap(pref)

% Get preferences:
if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1, 2); % Pre-allocate pass matrix.

dom = [-2 7];
map = bndfun.createMap(dom);

pass(1) = all(map.for([-1 1]) == dom) && all(map.inv(dom) == [-1 1]);

pass(2) = map.der(1) == diff(dom)/2;

end
