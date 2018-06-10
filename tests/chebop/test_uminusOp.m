function pass = test_uminusOp(~)
% TEST_UMINUSOP   Test CHEBOPs where unitary minus is involved. See #2191.
L = chebop(@(y) diff(y)-1,[0 1]);
L.lbc = 0;
y = L\0;
pass(1) = y(1) == 1;
%%
L.op = @(y) -diff(y)+1;
L.lbc = 0;
y = L\0;
pass(2) = y(1) == 1;
%%
L.op = @(y) -1*diff(y)+1;
y=L\0;
pass(3) = y(1) == 1;

%%
L.op = @(y) 1-diff(y);
y=L\0;
pass(4) = y(1) == 1;