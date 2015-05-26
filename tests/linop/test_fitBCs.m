function pass = test_fitBCs(pref)

if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

% Loop over discretization types:
type = {@chebcolloc1, @chebcolloc2, @ultraS};
err = Inf + zeros(length(type), 4);
for k = 1:3
    
    % Set disc type:
    pref.discretization = type{k};
    
    [Z, I, D, C, M] = linop.primitiveOperators();
    [zr, ev, su, dt] = linop.primitiveFunctionals();

    %% Scalar equation:
    bc = @(x, u) [sum(u) - 2 ; u(-1) + u(1) - 3];
    L = linop(D^2);
    L = addConstraint(L, su, -2);
    L = addConstraint(L, ev(-1) + ev(1), -3);
    u = fitBCs(L, pref);
    x = chebfun('x');
    err(k,1) = norm(bc(x, u{1}));
        
    %% Scalar equation (higher-order constraints):
    bc = @(x, u) [sum(u) - 2 ; u(-1) + u(1) - 3 ; feval(diff(u, 2), -0)];
    L = linop(D^3);
    L = addConstraint(L, su, -2);
    L = addConstraint(L, ev(-1) + ev(1), -3);
    L = addConstraint(L, ev(0)*D^2, 0);
    u = fitBCs(L, pref);
    x = chebfun('x');
    err(k,2) = norm(bc(x, u{1}));
    
    %% System:
    bc = @(x, u, v) [sum(u) - 1 ; u(0) - v(1) ; v(.5) - v(.7) ; sum(v) - 2];
    L = linop([D^2 I ; D^2 -I]);
    L = addConstraint(L, [su zr], -1);
    L = addConstraint(L, [ev(0) -ev(1)], 0);
    L = addConstraint(L, [zr ev(.5)-ev(.7)], 0);
    L = addConstraint(L, [zr su], -2);
    uu = fitBCs(L, pref);
    x = chebfun('x');
    err(k,3) = norm(bc(x, uu{1}, uu{2}));
    
    %% System with jumps:
    dom = [0 0.3 1];
    [Z, I, D, C] = linop.primitiveOperators(dom);
    [z, e, s] = linop.primitiveFunctionals(dom);
    j = functionalBlock.jumpAt(dom);
    A = linop( [ D^2, I; -D D^2+I ] );
    % Add constraints:
    A = addConstraint(A,[e(0) z],-1);
    A = addConstraint(A,[e(1) z],0);
    A = addConstraint(A,[z e(0)],0);
    A = addConstraint(A,[z e(1)],1);
    A = addContinuity(A,[j(0.3,1) z],2);
    A = addContinuity(A,[z j(0.3,1)],0);
    A = addContinuity(A,[j(0.3,0) z],0);
    A = addContinuity(A,[z j(0.3,0)],-.5);
    uu = fitBCs(A, pref);
    u1 = uu{1};
    u2 = uu{2};
    err(k, 4) = norm([u1(0)-1
                u1(1)
                u2(0)
                u2(1)+1
                jump(diff(u1), .3)+2
                jump(u1, .3)
                jump(u2, .3)-.5
                jump(diff(u2), .3)]);

end

pass = err < tol;
    
end
