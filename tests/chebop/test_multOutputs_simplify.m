function pass = test_multOutputs_simplify(~)
% Check that simplify gets called at the end of solving systems

%% Call with multiple outputs
% Linear chebop with two variables
L = chebop(@(t,u,v) [diff(u,2)+u; diff(v,2)+100*v],[0,10]);
L.lbc = @(u,v) [u-1; diff(u); v-1; diff(v)];
% Solve with two outputs
[u,v] = L\0;
% Their lengths should not be the same
pass(1) = length(u) ~= length(v);

%% Call with single output, should still be of same length
uv = L\0;
pass(2) = length(uv{1}) == length(uv{2});

%% Compare with previous outputs using deal
[u2, v2] = deal(uv);

pass(3) = length(u) == length(u2) && length(v) == length(v2);

end