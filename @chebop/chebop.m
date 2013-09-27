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
            
            numVars = nargin(N.op) - 1;
            
            if ( nargin < 2 )
                x = chebfun(@(x) x, N.domain);
            end
            if ( nargin < 3 )
                % Initialise a zero ADCHEBFUN:
                zeroFun = chebfun(0, N.domain);
                u = cell(numVars, 1);
                for k = 1:numVars
                    u{k} = zeroFun;
                end
                u = chebmatrix(u);
            end
            
            % Initialise the dependent variable:
            x = chebfun(@(x) x, N.domain);

            % Linearise
            [L, affine, isLinear] = linearise(N, x);
            
            %%
            % Solve:
            if ( isLinear )
                u = L\(rhs - affine);
            else
                if ( ~isempty(N.init) )
                    u = N.init;
                    [L, res] = linearise(N, x, u);
                    du = L\(rhs-res);
                else
                    du = L\rhs;
                end
                u = u - du;
                ub = u.blocks;
                res = N.op(x, ub{:}) - rhs;
                for newt = 1:10
                    % Linearise around current solution:
                    L = linearise(N, x, ub, flag);
                    % Solve the linearised system:
                    du = L\res;
                    % Append the Newton step:
                    u = u - du;
                    % Evaluate the residual
                    ub = u.blocks;
                    res = N.op(x, ub{:}) - rhs;
                                        
                    % Stop if well converged, or stagnated:
                    normStep(newt) = sum(cellfun(@norm, du.blocks));
                    if ( isa(res, 'chebmatrix') )
                        normRes(newt) = sum(cellfun(@norm, res.blocks));
                    else
                        normRes(newt) = norm(res);
                    end
                    if ( normStep(newt) < 1e-12 )
                        break
                    elseif (newt > 3 && normStep(newt) > 0.1*normStep(newt-3))
                        warning('CHEBFUN:bvpsc','Newton iteration stagnated.')
                        normStep.'
                        normRes.'
                        break
                    end
                    
                end
            end
            
        end
        
        function [L, res, isLinear] = linearise(N, x, u, flag)
            isLinear = true;
            numVars = nargin(N.op) - 1;
            
            if ( nargin < 2 )
                x = chebfun(@(x) x, N.domain);
            end
            if ( nargin < 3 )
                % Initialise a zero ADCHEBFUN:
                zeroFun = chebfun(0, N.domain);
                u = cell(numVars, 1);
                for k = 1:numVars
                    u{k} = zeroFun;
                end
%                 u0 = chebmatrix(u0);
            end

            if ( isa(u, 'chebmatrix') )
                u = u.blocks;
            end
            if ( isa(u, 'chebfun') )
                u = {u};
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
                if ( nargin == 4 )
                    lbcU.func = -lbcU.func;
                end
                L = bc(L, lbcU.jacobian, -lbcU.func); %#ok<CPROP>
                isLinear = isLinear & all(lbcU.isConstant);
            end

            if ( ~isempty(N.rbc) )
                rbcU = feval(N.rbc(u{:}), N.domain(end));
                if ( nargin == 4 )
                    rbcU.func = -rbcU.func;
                end
                L = bc(L, rbcU.jacobian, -rbcU.func); %#ok<CPROP>
                isLinear = isLinear & all(rbcU.isConstant);
            end

            if ( ~isempty(N.bc) )
                bcU = N.bc(x, u{:});
                if ( nargin == 4 )
                    bcU.func = -bcU.func;
                end
                L = bc(L, bcU.jacobian, -bcU.func); %#ok<CPROP>
                isLinear = isLinear & all(bcU.isConstant);
            end
            
        end
        
    end
end

