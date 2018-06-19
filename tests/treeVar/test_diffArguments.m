function pass = test_diffArguments(~)
% TEST_DIFFARGUMENTS   Test converting functions with more complicated
%                      arguments inside the diff method. These should all throw
%                      an error, see #2191.

%% Setup
dom = [-1, 4];
x = chebfun(@(x) x, dom);

try
    treeVar.toFirstOrder(@(u) diff(-u)+u, 0, dom);
catch ME
    pass(1) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:diff:diffArguments');
end

try
    treeVar.toFirstOrder(@(u) diff(2*u)+u, 0, dom);
catch ME
    pass(2) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:diff:diffArguments');
end

try
    treeVar.toFirstOrder(@(u) diff(2*u,2)+u, 0, dom);
catch ME
    pass(3) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:diff:diffArguments');
end

try
    treeVar.toFirstOrder(@(u) diff(x.*u)+u, 0, dom);
catch ME
    pass(4) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:diff:diffArguments');
end

try
    treeVar.toFirstOrder(@(u) diff(u+5,2)+u, 0, dom);
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:diff:diffArguments');
end