function pass = test_linjump

% Test jump condition in an ODE.
% TAD, 31 Jan 2014

tol = 1e-8;

solver = { @colloc2, @colloc1, @ultraS }; % FIXME
kind = [2 1 2];

for k = 1:length(solver)
    
    domain = [0 0.3 1];
    [Z,I,D,C] = linop.primitiveOperators(domain);
    [z,e,s] = linop.primitiveFunctionals(domain);
    j = functionalBlock.jumpAt(domain);

    % ODE:
    A = linop( [ D^2, I; -D D^2+I ] );

    % Add constraints:
    A = addConstraint(A,[e(0) z],-1);
    A = addConstraint(A,[e(1) z],0);
    A = addConstraint(A,[z e(0)],0);
    A = addConstraint(A,[z e(1)],1);
    A = addContinuity(A,[j(0.3,1) z],2);
    A = addContinuity(A,[z j(0.3,1)],0);
    A = addContinuity(A,[j(0.3,0) z],0);
    A = addContinuity(A,[z j(0.3,0)],0);

    x = chebfun('x',domain,'chebkind',kind(k));
    prefs = cheboppref;
    prefs.discretization = solver{k};
    u = linsolve(A,[x;0*x],prefs);

    %%
    % jumps
    J = functionalBlock.jump(0.3,domain,0);
    err(k,1) = J*u{1} + J*u{2};
    err(k,2) = J*(D*u{1}) - 2;

    %%
    % BCs
    err(k,3) = feval(u{1},0) + 1;
    err(k,4) = feval(u{2},1) - 1;

    %%
    % ODEs
    % [TODO]: The following lines use a hack. This needs to be changed:
    g = D^2*u{1}+u{2};
    g.funs{1} = g.funs{1}.funPart;
    g.funs{2} = g.funs{2}.funPart;
    err(k,5) = norm(g - x);
    % TODO: We would like to do this but deltafuns interfere:
    % err(k,5) = norm( D^2*u{1} + u{2} - x);
    err(k,6) = norm( -D*u{1} + D^2*u{2} + u{2} ); 
    
end

pass = err < tol;

end