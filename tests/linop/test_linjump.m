function pass = test_linjump

% Test jump condition in an ODE.
% TAD, 31 Jan 2014

tol = 1e-8;

solver = { @chebcolloc2, @chebcolloc1, @ultraS };
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
    A = addContinuity(A,[z j(0.3,0)],1);

    x = chebfun('x',domain,'chebkind',kind(k));
    prefs = cheboppref;
    prefs.discretization = solver{k};
    u = linsolve(A, [x ; 0*x], prefs);

    %%
    % jumps
    J = functionalBlock.jump(0.3,domain,0);
    err(k,1) = abs(J*u{2} - 1);
    err(k,2) = abs(J*(D*u{1}) - 2);

    %%
    % BCs
    err(k,3) = abs(feval(u{1},0) + 1);
    err(k,4) = abs(feval(u{2},1) - 1);

    %%
    % ODEs
    
    % TODO: We would like to do this but deltafuns interfere:
    % err(k,6) = norm( D^2*u{1} + u{2} - x);
    
    
    res = D^2*u{1} + u{2} - x;
    if ( isa(res.funs{1}, 'deltafun') )
        res.funs{1} = res.funs{1}.funPart;
    end
    if ( isa(res.funs{2}, 'deltafun') )
        res.funs{2} = res.funs{2}.funPart;
    end    
    err(k,5) = norm( res);

    
    res = -D*u{1} + D^2*u{2} + u{2};
    if ( isa(res.funs{1}, 'deltafun') )
        res.funs{1} = res.funs{1}.funPart;
    end
    if ( isa(res.funs{2}, 'deltafun') )
        res.funs{2} = res.funs{2}.funPart;
    end    
    err(k,6) = norm( res ); 
    
end

pass = err < tol;

end
