%CHEBOP  Construct an operator on chebfuns.

classdef (InferiorClasses = {?double}) chebop
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        domain = [];    % Domain of the operator
        op = [];        % The operator
        lbc = [];       % Left boundary condition(2)
        rbc = [];       % Right boundary condition(s)
        bc = [];        % Other/internal/mixed boundary conditions
        init = [];      % Initial guess of solution
    end
    
    methods
        
        function N = chebop(op, dom)
            if ( nargin == 0 )
                return
            end
            if ( nargin < 2 )
                dom = chebfun.pref('domain');
            end
            
            N.op = op;
            N.domain = dom;
            
        end
        
        function u = mldivide(N, rhs)
            
            %%
            % Initialise an ADCHEBFUN:
            numVars = nargin(N.op) - 1;
            zeroFun = chebfun(0, N.domain);
            u0 = cell(numVars, 1);
            for k = 1:numVars
                u0{k} = zeroFun;
            end

            %%
            % Initialise the dependent variable:
            x = chebfun(@(x) x, N.domain);

            [L, res, isLinear] = linearise(N, u0, x);
            
            %%
            % Solve:
            u = L\(rhs-res);
            
            if ( ~isLinear )
                ub = u.blocks;
                res = N.op(x, ub{:}) - rhs;
                for newt = 1:10
                    % Linearise around current solution:
                    L = linearise(N, ub, x);
                    % Solve the linearised system:
                    du = L\res;
                    % Append the Newton step:
                    u = u - du;
                    % Evaluate the residual
                    ub = u.blocks;
                    res = N.op(x, ub{:}) - rhs;
                                        
                    % Stop if well converged, or stagnated:
                    normStep(newt) = sum(cellfun(@norm, du.blocks));
                    normRes(newt) = sum(cellfun(@norm, res.blocks));
                    if ( normStep(newt) < 1e-12 )
                        break
                    elseif (newt > 3 && normStep(newt) > 0.1*normStep(newt-3))
                        warning('CHEBFUN:bvpsc','Newton iteration stagnated.')
                        break
                    end
                    
                end
            end
            
        end
        
        function [L, res, isLinear] = linearise(N, u, x)
            isLinear = true;
            numVars = nargin(N.op) - 1;
            
            if ( isa(u, 'chebmatrix') )
                u = u.blocks;
            end
            
            for k = 1:numVars
                u{k} = seed(adchebfun(u{k}), k, numVars);
            end
            
            
            %%
            % Evaluate the operators to get a linearisation:
            w = N.op(x, u{:});
            L = vertcat(w.jacobian);
            res = vertcat(w.func);
            isLinear = isLinear & all(w.isConstant);

            %%
            % Add BCs
            if ( ~isempty(N.lbc) )
                lbcU = feval(N.lbc(u{:}), N.domain(1));
                L = bc(L, lbcU.jacobian, -lbcU.func); %#ok<CPROP>
                isLinear = isLinear & all(lbcU.isConstant);
            end

            if ( ~isempty(N.rbc) )
                rbcU = feval(N.rbc(u{:}), N.domain(end));
                L = bc(L, rbcU.jacobian, -rbcU.func); %#ok<CPROP>
                isLinear = isLinear & all(rbcU.isConstant);
            end

            if ( ~isempty(N.bc) )
                bcU = N.bc(x, u{:});
                L = bc(L, bcU.jacobian, -bcU.func); %#ok<CPROP>
                isLinear = isLinear & all(bcU.isConstant);
            end
            
        end
        
    end
end

