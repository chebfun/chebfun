function pass = test_bivariate(~)
%TEST_BIVARIATE   Test various univariate treeVar methods
%   This methods generates a number of anonymous functions that we convert to
%   first order format. We then evaluate them, and check that they match the
%   expected results.

%% Setup
% Initialize a TREEVAR variable.
dom = [0 2];

% Used for setting up anonymous functions we convert to first order form.
rhs = 3.2;
alp = 2;

% Arguments we eventually evaluate the functions at:
tArg = 1;
uArg = [.4 .2];
sArg = 2;

% List of methods that we want to test:
testMethods = {@minus, @plus, @power, @rdivide, @times, @mrdivide, @mtimes};

% Store comparison errors:
errors = zeros(length(testMethods), 1);
for mCounter = 1:length(testMethods)-2
    % Current method we're testing:
    method = testMethods{mCounter};
    
    % Construct an anonymous function that calls the current method of
    % interest:
    myFun = @(u) diff(u, 2) + diff(u) +  alp*method(u, diff(u));
    
    % Convert MYFUN to first order system:
    anonFun = treeVar.toFirstOrder(myFun, rhs, dom);
    
    % The correct first order reformulation of MYFUN:
    correctFun = @(t,u) [u(2); rhs - u(2) - alp*method(u(1), u(2))];
    
    % Compare the result of evaluating the automatically converted ANONFUN and
    % the manually constructed CORRECTFUN. The difference in the outputs
    % should be very small (if any):
    errors(mCounter) = norm(anonFun(tArg, uArg) - correctFun(tArg, uArg));
end

for mCounter = length(testMethods)-1:length(testMethods)
    % Current method we're testing:
    method = testMethods{mCounter};
    
    % Construct an anonymous function that calls the current method of
    % interest:
    myFun = @(u) diff(u, 2) + diff(u) +  alp*method(u, sArg);
    
    % Convert MYFUN to first order system:
    anonFun = treeVar.toFirstOrder(myFun, rhs, dom);
    
    % The correct first order reformulation of MYFUN:
    correctFun = @(t,u) [u(2); rhs - u(2) - alp*method(u(1), sArg)];
    
    % Compare the result of evaluating the automatically converted ANONFUN and
    % the manually constructed CORRECTFUN. The difference in the outputs
    % should be very small (if any):
    errors(mCounter) = norm(anonFun(tArg, uArg) - correctFun(tArg, uArg));
end

pass = errors < 10*eps;
