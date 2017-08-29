function pass = test_ivp_chebmatrix_syntax
% Test CHEBMATRIX syntax for IVPs, using full Brusselator IVP:

%% Setup
dom = [0 5];
tol = 1e-14;

%% With CHEBOP, multiple variable syntax for reference
N = chebop(@(t,u,v,w) [
    diff(u)-1-u^2*v + (w+1)*u;
    diff(v) - u*w+u^2*v;
    diff(w) - 1.5 + w*u], dom);
%% Test with left and right BCs:
N.lbc = [1;3;4];
uL = N\0;
N.lbc = [];
pass(1) = norm(feval(uL,dom(1)) - [1;3;4]) < tol;
%%
N.lbc = [];
N.rbc = [.7; 0.9; 1.2];
uR = N\0;
pass(2) = norm(feval(uR,dom(end)) - [.7; 0.9; 1.2]) < tol;
%% CHEBMATRIX syntax

% Define a new CHEBOP
M = chebop(@(t,u) [
    diff(u{1})-1-u{1}^2*u{2} + (u{3}+1)*u{1};
    diff(u{2}) - u{1}*u{3}+u{1}^2*u{2};
    diff(u{3}) - 1.5 + u{3}*u{1}], dom);

%% CHEBMATRIX left boundary conditions
M.lbc = @(u)[u{1}-1; u{2}-3; u{3}-4];
vL1 = M\0;
% Results should be identical to above
pass(3) = norm(vL1 - uL) == 0;

%% Vector LBC
M.lbc = [1; 3; 4];
vL2 = M\0;
pass(4) = norm(vL2 - uL) == 0;

%% CHEBMATRIX right boundary conditions
M.lbc = [];
M.rbc = @(u)[u{1}-.7; u{2}-0.9; u{3}-1.2];
vR1 = M\0;
% Results should be identical to above
pass(5) = norm(vR1 - uR) == 0;

%% Vector RBC
M.rbc = [.7; .9; 1.2];
vR2 = M\0;
pass(6) = norm(vR2 - uR) == 0;

end

