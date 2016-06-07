function pass = test_linearizationDimensions
% TEST_LINEARIZATIONDIMENSIONS    Test that we get the correct dimensions when
% linearizing parameter dependent problems. CF #1176.

%% Example 1
L = chebop(@(x,u,v) u + v);
LL = linop(L); LLB = LL.blocks;
% Should expect both blocks to be OPERARTORBLOCK, since no variable gets
% differentiated or integrated.
pass(1) = isa(LLB{1}, 'operatorBlock') && isa(LLB{2}, 'operatorBlock');

%% Example 2:
L = chebop(@(x,u) u{1} + u{2});
LL = linop(L); LLB = LL.blocks;
% Should expect both blocks to be OPERARTORBLOCK, since no variable gets
% differentiated or integrated.
pass(2) = isa(LLB{1}, 'operatorBlock') && isa(LLB{2}, 'operatorBlock');

%% Example 3:
L = chebop(@(x,u,v) diff(u) + v);
LL = linop(L); LLB = LL.blocks;
% Should expect both the first block to be an OPERATORBLOCK, the second to be a
% chebfun.
pass(3) = isa(LLB{1}, 'operatorBlock') && isa(LLB{2}, 'chebfun'); 
%% Example 4:
L = chebop(@(x,u,v) u + diff(v));
LL = linop(L); LLB = LL.blocks;
% Here, we should expect both the first block to be an OPERATORBLOCK, the second
% to be a chebfun.
pass(4) = isa(LLB{1}, 'chebfun') && isa(LLB{2}, 'operatorBlock'); 
%% Example 5:
L = chebop(@(x,u,v) diff(u) + diff(v));
LL = linop(L); LLB = LL.blocks;
% Here, we should expect both blocks to be an ç.
pass(5) = isa(LLB{1}, 'operatorBlock') && isa(LLB{2}, 'operatorBlock');
%% Example 6:
L = chebop(@(x,u,v) [u + diff(v); diff(u) + v]);
LL = linop(L); LLB = LL.blocks;
% Here, we should expect both the first block to be an operatorBlock, the second
% to be a chebfun.
pass(6) = isa(LLB{1,1}, 'operatorBlock') && isa(LLB{2,1}, 'operatorBlock')  && ...
    isa(LLB{1,2}, 'operatorBlock') && isa(LLB{2,2}, 'operatorBlock');

%% Example 7 -- more than two variables
L = chebop(@(x,u,v,w) u + diff(v) + w);
LL = linop(L); LLB = LL.blocks;
% Here, we should expect both the first block to be an operatorBlock, the second
% to be a chebfun.
pass(7) = isa(LLB{1}, 'chebfun') && isa(LLB{2}, 'operatorBlock') && ...
    isa(LLB{3}, 'chebfun');

%% Example 8 -- more than two variables
L = chebop(@(x,u,v,w) [u + diff(v) + w; u + v + diff(w)]);
LL = linop(L); LLB = LL.blocks;
% Here, we should expect both the first block to be an operatorBlock, the second
% to be a chebfun.
pass(8) = isa(LLB{1,1}, 'chebfun') && isa(LLB{2,1}, 'chebfun') &&  ...
    isa(LLB{1,2}, 'operatorBlock') && isa(LLB{2,2}, 'operatorBlock') && ...
    isa(LLB{1,3}, 'operatorBlock') && isa(LLB{2,3}, 'operatorBlock');

end