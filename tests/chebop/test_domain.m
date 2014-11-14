function pass = test_domain(pref)
%TEST_DOMAIN    Ensure that we can both pass row and column vectors to CHEBOP

%% Setup
if ( nargin == 0)
    pref = cheboppref();
end

dom = [0,2];

%% Row vector passed to constructor
N = chebop(@(u) diff(u, 2) + u, dom);
N.bc = 0;
u1 = N\1;

%% Column vector passed to constructor
N = chebop(@(u) diff(u, 2) + u, dom');
N.bc = 0;
u2 = N\1;

%% Row vector set after construction
N = chebop(@(u) diff(u, 2) + u);
N.domain = dom;
N.bc = 0;
u3 = N\1;

%% Column vector set after construction
N = chebop(@(u) diff(u, 2) + u);
N.domain = dom';
N.bc = 0;
u4 = N\1;

%% Check that all solutions are the same
pass = ( norm(u1 - u2) + norm(u3 - u4) + norm(u1 - u4) ) == 0;

end
