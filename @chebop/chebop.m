classdef (InferiorClasses = {?double}) chebop
%CHEBOP  CHEBOP class for representing operators on functions defined on [a,b].
   
    properties ( GetAccess = 'public', SetAccess = 'public' )
        domain = [];    % Domain of the operator
        op = [];        % The operator
        lbc = [];       % Left boundary condition(s)
        rbc = [];       % Right boundary condition(s)
        bc = [];        % Other/internal/mixed boundary conditions
        init = [];      % Initial guess of a solution
        % Default discretization for linear problems
        discretizationType = @colloc2; 
    end
    
    methods
        
        function N = chebop(op, dom)
            if ( nargin == 0 )
                return
            end
            if ( nargin < 2 )
                % Need to access chebfunpref to create an operator on the
                % default domain if none is passed.
                p = chebpref();
                dom = p.domain;
            end
            
            N.op = op;
            N.domain = dom;
            
        end
        
        function u = mldivide(N, rhs)
            pref = cheboppref;
            u = solvebvp(N, rhs, pref);
        end
            
        [L, res, isLinear] = linearise(N, x, u, flag);
        
        
    end
        
    methods (Static = true) % These should be private methods as well
        
        [displayFig, displayTimer] = displayInfoInit(u,pref);
        
        displayInfoIter(u, delta, iterNo, normdu, cFactor, lendu, ...
            lambda, lenu, displayFig, displayTimer, pref);
        
        displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, ...
            displayTimer, pref)
        
    end
    
end

