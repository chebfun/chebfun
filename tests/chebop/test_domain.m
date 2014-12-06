function pass = test_domain(~)
%TEST_DOMAIN    Ensure that we can both pass row and column vectors to CHEBOP

%% Setup
dom = [0,2];

%% Column vector passed to constructor -- should give an error
try
    chebop(@(u) diff(u, 2) + u, dom');
catch ME
   pass(1) = strcmp(ME.identifier, 'CHEBOP:CHEBOP:domain');
end

%% Column vector set after construction -- should give an error
N = chebop(@(u) diff(u, 2) + u);
try
    N.domain = dom';
catch ME
   pass(2) = strcmp(ME.identifier, 'CHEBOP:SET:domain');
end

%% Row vector passed to constructor -- should be OK
N = chebop(@(u) diff(u, 2) + u, dom);
N.bc = 0;
u1 = N\1;

%% Row vector set after construction -- should be OK
N = chebop(@(u) diff(u, 2) + u);
N.domain = dom;
N.bc = 0;
u2 = N\1;

%% Check that both allowed approaches give the same solution
pass(3) = norm(u1 - u2) == 0;

end
