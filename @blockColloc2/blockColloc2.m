classdef blockColloc2 < blockDiscretization
    properties 
        size = [];  % arbitrary, but fixed in any one instance
        domain = [-1 1];
    end
    
    methods
        function A = blockColloc2(varargin)
            % Collocation matrix on 2nd kind points.
            
            % COLLOC2(DIM,DOMAIN) returns a dummy object that will propagate the
            % dimension size DIM and function domain DOM throughout the delayed
            % evaluation stack.
            %
            % COLLOC2(A,DIM) realizes the linop A (which knows its domain) at
            % dimension DIM.
            % 
            % COLLOC2([]) returns a dummy object that gives access to static
            % methods. 
            
            if ( nargin > 1 )
                if isa(varargin{1}, 'linBlock')
                    L = varargin{1};
                    A.size = varargin{2};
                    A.domain = L.domain;
                    A = L.delayFun( A );
                else
                    A.size = varargin{1};
%                     validateattributes(varargin{2},{'numeric'},{'increasing','finite'});
                    A.domain = varargin{2};
                end
            end
        end
                
        % Building blocks.
        
        D = diff(A,m) 
        C = cumsum(A,m)
        
        function I = eye(A)
            n = dim(A);
            I = eye(sum(n));
        end
        
        function Z = zeros(A)
            n = dim(A);
           Z = zeros(sum(n));
        end
        
        function Z = zero(A)
            n = dim(A);
            Z = zeros(1,sum(n));
        end
        
        F = diag(A,f)
        
        % Required operators.
        function C = mtimes(A,B)
            C = A*B;
        end
        
        function C = plus(A,B)
            C = A+B;
        end
        
        function B = uminus(A)
            B = -A;
        end
        
        % Required functionals.
        
        S = sum(A)
        
        E = feval(A,location,direction)
        
        function F = inner(A,f)
            d = chebmatrix.mergeDomains({A,f});
            [x,w] = points(dim(A),d,2);
            F = w.*f(x);
        end
        
    end
    
    methods (Static)
        % Additional methods
        
        B = resize(A,m,n,domain)
        
        [isDone,epsLevel] = testConvergence(v)
        
        function fx = discretizeFunction(f,dim,dom)
            if ( nargin < 3 )
                dom = f.domain;
            end
            x = blockColloc2.points(dim,dom);
            fx = f(x);
        end
        
        [x,w] = points(n,d)
        
                
        function L = discretize(A, dim, dom)
            L = A.delayFun( blockColloc2(dim, dom) );
        end
        
        function f = makeChebfun(u, dom)
            f = chebfun(u, dom);
        end
        
    end
    
    methods (Access=private)
        
        D = diffmat(N,k)
       
        Q = cumsummat(N)
                
    end
end