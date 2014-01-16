% Test file for fun/isempty.m

function pass = test_isempty(varargin)

% Set a domain for BNDFUN.
dom = [-2 7];

%%
% Test an empty BNDFUN:
f = bndfun();
pass(1) = isempty(f);
    
%%
% Test an non-empty BNDFUN:
f = bndfun(@sin, dom); 
pass(2) = ~isempty(f);

%%
% Test an non-empty array-valued BNDFUN:
f = bndfun(@(x) [sin(x), cos(x)], dom);
pass(3) = ~isempty(f);
    
%%
% Test an array of BNDFUN objects:
f = [ bndfun(@sin, dom), bndfun(@sin, dom) ];
pass(4) = ~isempty(f);

%% 
% [TODO]: Run a few tests for UNBNDFUN.
end