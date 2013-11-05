classdef blockFunction
    
    properties (Access=public)
        % This property is assigned the callable function that does the
        % correct operation when called on chebfuns.
        func = [];
    end
    
    methods

        function A = blockFunction(f)
            A.func = f;  
        end
        
        % ** Implementation of abstract methods **
        
        % Required operators.
        
         
        % ** Additional operations **
        
        function C = mtimes(A,B)
            C = blockFunction( @(z) A.func(B.func(z)) );
        end
        
        function C = plus(A,B)
            C = blockFunction( @(z) A.func(z) + B.func(z) );
        end
        
        function C = uminus(A)
            C = blockFunction( @(z) -A.func(z) );
        end
        
    end
    
    methods (Static)
        function D = diff(domain,m)
            D = blockFunction( @(z) diff(z,m) );
        end
        
        function C = cumsum(domain,m)
            C = blockFunction( @(u) cumsum(u,m) );
        end
        
        function I = eye(domain)
            I = blockFunction( @(z) z );
        end
        
        function Z = zeros(domain)
            Z = blockFunction( @(z) chebfun(0,domain) );
        end
        
        function F = mult(f)
            F = blockFunction( @(z) times(f,z) );
        end    
                
        function S = sum(domain)
            S = blockFunction( @(z) sum(z) );
        end
        
        function E = feval(domain,location,direction)
            if (direction < 0)
                E = blockFunction( @(u) feval(u,location,'left') );
            elseif (direction > 0)
                E = blockFunction( @(u) feval(u,location,'right') );
            else
                E = blockFunction( @(u) feval(u,location) );
            end
        end
        
        function F = inner(f)
            F = blockFunction( @(z) mtimes(f',z) );
        end


    end
    
end
        