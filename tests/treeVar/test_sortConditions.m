function pass = test_sortConditions(~)
%TEST_SORTCONDITIONS   Test that we sort IC/FCs correctly.

%% Setup
dom = [.4 2];

%% First order scalar IVP
icFun = @(u) u;
correct = 1;
idx = treeVar.sortConditions(icFun, dom, 1);
pass(1) = ( idx == correct );

%% Second order scalar IVP, increasing order
icFun = @(u) [u; diff(u)];
correct = [1 2];
idx = treeVar.sortConditions(icFun, dom, 2);
pass(2) = all( idx == correct );

%% Second order scalar IVP, decreasing order
icFun = @(u) [diff(u); u];
correct = [2 1];
idx = treeVar.sortConditions(icFun, dom, 2);
pass(3) = all( idx == correct );

%% Sixth order scalar IVP, increasing order
icFun = @(u) [diff(u, 0); diff(u, 1); diff(u,2); diff(u,3); diff(u,4); diff(u,5)];
correct = 1:6;
idx = treeVar.sortConditions(icFun, dom, 6);
pass(4) = all( idx == correct );

%% Sixth order scalar IVP, arbitrary order
icFun = @(u) [diff(u, 2); diff(u, 4); diff(u,5); diff(u,0); diff(u,1); diff(u,3)];
correct = [4 5 1 6 2 3];
idx = treeVar.sortConditions(icFun, dom, 6);
pass(5) = all( idx == correct );

%% First order coupled IVP, two variables, increasing order
icFun = @(u,v) [u; v];
correct = [1 2];
idx = treeVar.sortConditions(icFun, dom, [1 1]);
pass(6) = all( idx == correct );

%% First order coupled IVP, two variables, decreasing order
icFun = @(u,v) [v; u];
correct = [2 1];
idx = treeVar.sortConditions(icFun, dom, [1 1]);
pass(7) = all( idx == correct );

%% Second order coupled IVP, two variables, increasing order
icFun = @(u,v) [u; diff(u); v; diff(v)];
correct = [1 2 3 4];
idx = treeVar.sortConditions(icFun, dom, [2 2]);
pass(8) = all( idx == correct );

%% Second order coupled IVP, two variables, arbitrary order
icFun = @(u,v) [u; diff(v); diff(u); v];
correct = [1 3 4 2];
idx = treeVar.sortConditions(icFun, dom, [2 2]);
pass(9) = all( idx == correct );

%% Fourth order coupled IVP, four variables, arbitrary order
icFun = @(u, v, w, y) [diff(u, 3); diff(w,3); w; diff(u,2); ...
    diff(v, 3); diff(v, 2); diff(w); y; ...
    diff(u); diff(y, 2); diff(v); v; ...
    diff(y, 1); diff(y, 3); u; diff(w, 2)];
correct = [15 9 4 1 12 11 6 5 3 7 16 2 8 13 10 14];
idx = treeVar.sortConditions(icFun, dom, [4 4 4 4]);
pass(10) = all( idx == correct );

%% Unsupported format, multiplying unknown function
icFun = @(x,u) 5*u -1;
try
    treeVar.sortConditions(icFun, dom, 1);
catch ME
    % The highest order derivatives of u and v appear in the same line -- this
    % should give us an error.
    errorPass(1) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:unsupportedCondition');
end

%% Unsupported format, multiplying derivative of unknown function
icFun = @(x,u) 5*diff(u) -1;
try
    treeVar.sortConditions(icFun, dom, 1);
catch ME
    % The highest order derivatives of u and v appear in the same line -- this
    % should give us an error.
    errorPass(2) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:unsupportedCondition');
end

%% Unsupported format, unknown function appears twice
icFun = @(x,u) u + diff(u);
try
    treeVar.sortConditions(icFun, dom, 2);
catch ME
    % The highest order derivatives of u and v appear in the same line -- this
    % should give us an error.
    errorPass(3) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:unsupportedCondition');
end

%% Unsupported format, unknown function appears twice, system
icFun = @(x,u,v) u + diff(u);
try
    treeVar.sortConditions(icFun, dom);
catch ME
    % The highest order derivatives of u and v appear in the same line -- this
    % should give us an error.
    errorPass(4) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:unsupportedCondition');
end


%% Unsupported format, coupled conditions
icFun = @(x,u,v) u + diff(v);
try
    treeVar.sortConditions(icFun, dom, [1 2]);
catch ME
    % The highest order derivatives of u and v appear in the same line -- this
    % should give us an error.
    errorPass(5) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:nonSeparated');
end

%% Unsupported format, two conditions on same variable
icFun = @(u) [u - 1; u - 2];
try
    treeVar.sortConditions(icFun, dom, 2);
catch ME
    errorPass(6) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:multipleConditionsSameVariable');
end
%% Unsupported format, specify values of a too high order derivative 
icFun = @(y) [y-1; diff(y,2)];
try
    treeVar.sortConditions(icFun, dom, 2);
catch ME
    errorPass(7) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:tooHighOrderCondition');
end
%% Unsupported format, need to specify values of solution and all appr. derivs. 
icFun = @(y) [y-1; diff(y,2)];
try
    treeVar.sortConditions(icFun, dom, 3);
catch ME
    errorPass(8) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:missingConditions');
end

%% Unsupported format, not many enough derivative values specified:
icFun = @(y) [y-1; diff(y)];
try
    treeVar.sortConditions(icFun, dom, 3);
catch ME
    errorPass(9) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:missingConditions');
end

%% Unsupported format, two conditions on same variable
icFun = @(u,v) [u - 1; v - 2; diff(u) - 2; v - 3];
try
    treeVar.sortConditions(icFun, dom, [2 1]);
catch ME
    errorPass(10) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:multipleConditionsSameVariable');
end
%% Unsupported format, specify values of a too high order derivative 
icFun = @(u,v) [u-1; diff(v,2); diff(u); v-3];
try
    treeVar.sortConditions(icFun, dom, [2 2]);
catch ME
    errorPass(11) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:tooHighOrderCondition');
end
%% Unsupported format, need to specify values of solution and all appr. derivs. 
icFun = @(u,v) [u-1; diff(v,2); diff(v); diff(u,2)];
try
    treeVar.sortConditions(icFun, dom, [3 3]);
catch ME
    errorPass(12) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:missingConditions');
end

%% Unsupported format, not many enough derivative values specified:
icFun = @(u,v) [u-1; diff(u); diff(v); v-3];
try
    treeVar.sortConditions(icFun, dom, [3 3]);
catch ME
    errorPass(13) = strcmp(ME.identifier, ...
        'CHEBFUN:TREEVAR:sortConditions:missingConditions');
end

%% Combine the information
pass = [pass, errorPass];
end
