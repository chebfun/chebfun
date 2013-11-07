classdef blockColloc2 < blockDiscretization
    
    methods
        function A = blockColloc2(varargin)
            % Collocation matrix on 2nd kind points.
            
            % COLLOC2(DIM, DOMAIN) returns a dummy object that will propagate the
            % dimension size DIM and function domain DOM throughout the delayed
            % evaluation stack.
            %
            % COLLOC2(A, DIM) realizes the linop A (which knows its domain) at
            % dimension DIM.
            %
            % COLLOC2([]) returns a dummy object that gives access to static
            % methods.
            
            if ( nargin > 1 )
%                if isa(varargin{1}, 'linBlock')
%                     L = varargin{1};
%                     A.size = varargin{2};
%                     A.domain = L.domain;
%                     A = L.stack( A );
%                 else
                    A.domain = varargin{2};
                    A.dimension = varargin{1};
%                end
            end
        end
        
        % Building blocks.
        
        D = diff(A, m)
        C = cumsum(A, m)
        
        function I = eye(A)
            n = A.dimension;
            I = eye(sum(n));
        end
        
        function Z = zeros(A)
            n = A.dimension;
            Z = zeros(sum(n));
        end
        
        function Z = zero(A)
            n = A.dimension;
            Z = zeros(1, sum(n));
        end
        
        F = mult(A, f)
        
        % Required operators.
        function C = mtimes(A, B)
            C = A*B;
        end
        
        function C = plus(A, B)
            C = A+B;
        end
        
        function B = uminus(A)
            B = -A;
        end
        
        % Required functionals.
        
        S = sum(A)
        
        E = feval(A, location, direction)
        
        function F = inner(A, f)
            d = chebmatrix.mergeDomains({A, f});
            [x, w] = points(A);
            F = w.*f(x);
        end
        
    end
    
    methods
        % Additional methods
        
        %B = resize(A, m, n, domain, difforder)
        
        %[isDone, epsLevel] = testConvergence(v)
        
        function fx = toValues(disc,f)
            n = disc.dimension;

            x = points(disc);
            fx = f(x);

            % Evaluate left- and right-sided limits at breaks:
            csn = [0, cumsum(n)];
            dxloc = csn(2:end-1);
            fx(dxloc) = feval(f, x(dxloc), 'left');
            fx(dxloc+1) = feval(f, x(dxloc), 'right');
        end
               
        function L = matrix(disc,A)
            validateParameters(disc);
            if isa(A,'linBlock')
                L = A.stack( disc );
            elseif isa(A,'chebfun')
                L = disc.toValues(A);
                if ( A.isTransposed )
                    L = L.';
                end
            elseif isnumeric(A)
                L = A;
            else
                error('Unrecognized block type.')
            end
        end
        
        function f = toFunction(disc,values)
            f = chebfun(values, disc.domain);
        end
        
    end
    
end