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
        
        function D = diff(A,m)
            D = op( @(z) diff(z,m) );
        end
        
        function C = cumsum(A,m)
            C = op( @(u) cumsum(u,m) );
        end
        
        function I = eye(A)
            I = op( @(z) z );
        end
        
        function Z = zeros(A)
            Z = op( @(z) chebfun(0,domain) );
        end
        
        function F = diag(A,f)
            F = op( @(z) times(f,z) );
        end    
        
        % Required functionals.
        
        function S = sum(A)
            S = op( @(z) sum(z) );
        end
        
        function E = evalAt(A,location,direction)
            if (direction < 0)
                E = op( @(u) feval(u,location,'left') );
            elseif (direction > 0)
                E = op( @(u) feval(u,location,'right') );
            else
                E = op( @(u) feval(u,location) );
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
        