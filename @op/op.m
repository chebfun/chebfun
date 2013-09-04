classdef op 
    
    properties (Access=public)
        % This property is assigned the callable function that does the
        % correct operation when called on chebfuns.
        func = [];
    end
    
    methods

        function A = op(f)
            A.func = f;  
        end
        
        % ** Implementation of abstract methods **
        
        % Required operators.
        
        function D = diff(A,domain)
            D = op( @(z) diff(z) );
        end
        
        function C = cumsum(A,domain)
            C = op( @(u) cumsum(u) );
        end
        
        function I = eye(A,domain)
            I = op( @(z) z );
        end
        
        function Z = zeros(A,domain)
            Z = op( @(z) chebfun(0,domain) );
        end
        
        function F = diag(A,f)
            F = op( @(z) times(f,z) );
        end    
        
        % Required functionals.
        
        function S = sum(A,domain)
            S = op( @(z) sum(z) );
        end
        
        function E = evalAt(A,domain,loc)
            if strcmp(loc,'left')
                E = op( @(u) feval(u,domain(1)) );
            elseif strcmp(loc,'right')
                E = op( @(u) feval(u,domain(2)) );
            elseif isnumneric(loc)
                E = op( @(u) feval(u,loc) );
            end
        end
        
        function F = inner(A,f)
            F = op( @(z) mtimes(f',z) );
        end

        
        % ** Additional operations **
        
        function C = mtimes(A,B)
            C = op( @(z) A.func(B.func(z)) );
        end
        
        function C = plus(A,B)
            C = op( @(z) A.func(z) + B.func(z) );
        end
        
        function C = uminus(A)
            C = op( @(z) -A.func(z) );
        end
        
    end
    
end
        
