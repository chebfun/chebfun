function u = mldivide(N, rhs)

numVars = nargin(N.op) - 1;

% Initialise a zero ADCHEBFUN:
zeroFun = chebfun(0, N.domain);
u0 = cell(numVars, 1);
for k = 1:numVars
    u0{k} = zeroFun;
end
u0 = chebmatrix(u0);

% Initialise the dependent variable:
x = chebfun(@(x) x, N.domain);

% Linearise
[L, affine, isLinear] = linearise(N, x, u0);

    function out = mynorm(f)
        if ( isa(f, 'chebmatrix') )
            out = max(cellfun(@(u) get(u, 'vscale'), f.blocks));
        else
            out = get(f, 'vscale');
        end
    end

%%
% Solve:
if ( all(isLinear) )
    % Linear solve:
    u = L\(rhs - affine);
    
else
    
    % Newton solve:
    
    if ( ~isempty(N.init) )
        u = N.init;
        L = linearise(N, x, u);
        res = N.op(x, u{:}) - rhs;
        du = L\res;
        %                 else
        %                     u = makeGuess(L, x);
        %                     L = linearise(N, x, u);
        %                     res = N.op(x, u{:}) - rhs;
        %                     du = L\res;
    else
        u = u0;
        du = L\(rhs - affine);
    end
    
    u = u - du;
    ub = u.blocks;
    res = N.op(x, ub{:}) - rhs;
    
    fprintf('Nonlinear equation. using Newton''s method:\n')
    fprintf('step\t normUpdate\t\t  normRes\n')
    normUpdate(1,1) = mynorm(du);
    normRes(1,1) = mynorm(res);
    fprintf(' %2.2d\t%16.16f\t %16.16f\n', 1, normUpdate(1,1), normRes(1,1))
    newt = 1;
    while ( normRes(end) > 1e-12)
        newt = newt +1;
        % Linearise around current solution:
        L = linearise(N, x, ub, []); % flag to negate contraint RHSs.
        % Solve the linearised system:
        du = L\res;
        % Append the Newton step:
        u = u - du;
        % Evaluate the residual
        ub = u.blocks;
        res = N.op(x, ub{:}) - rhs;
        
        % Stop if well converged, or stagnated:
        normUpdate(newt,1) = mynorm(du);
        normRes(newt,1) = mynorm(res);
        fprintf(' %2.2d\t%16.16f\t %16.16f\n', newt, ...
            normUpdate(newt,1), normRes(newt,1))
        if ( normUpdate(newt,1) < 1e-12 )
            break
            %                     elseif (newt > 3 && normUpdate(newt) > 0.1*normUpdate(newt-3))
            %                         warning('CHEBFUN:bvpsc','Newton iteration stagnated.')
            %                         break
        elseif (newt > 10)
            warning('CHEBFUN:bvpsc','Newton iteration failed.')
            break
        end
        
        
    end
end

end