classdef operatorBlock < linBlock
    
    properties
    end
    
    methods
        function A = operatorBlock(domain)
            A = A@linBlock(domain);
        end
        
        function varargout = size(A, dim)
            % S = SIZE(A)
            % [M, N] = SIZE(A)
            % P = SIZE(A, K)
            %
            m = [Inf Inf];

            if ( nargin > 1 )
                varargout = {m(dim)};
            elseif ( nargout <= 1 )
                varargout = {m};
            else
                varargout = {m(1) m(2)};
            end
        end

        function C = uminus(A)
            C = operatorBlock(A.domain);
            C.stack = @(z) -A.stack(z);
            C.func = -A.func;
            C.coeff = -A.coeff;
            C.diffOrder = A.diffOrder; 
        end
        

        function C = mtimes(A, B)
            % A*B
            % If A, B both linops, or one is linop and one scalar, the
            % result is the composed linop.
            %
            % A*f
            % If A is a linop and f is a chebfun, return the chebfun
            % resulting from application of A to f.
            %
            % A*u
            % If A is a linop and u is a matrix, return the matrix
            % resulting from application of the discretization of A to u.
            %
            
            % No error checking here.
            % Which case?
            if ( isa(B, 'chebfun') )
                C = A.functionForm;
                C = C(B);
            elseif ( isnumeric(B) )
                N = size(B, 1);    % discretization size
                L = matrix(A, N);
                C = L*B;
            else
                % A scalar is converted into a constant chebfun, which is then
                % diagnified. 
                if ( isnumeric(A) )
                    A = linBlock.mult( chebfun(A, B.domain) );
                elseif ( isnumeric(B) )
                    B = linBlock.mult( chebfun(B, A.domain) );
                end
                
                C = operatorBlock(A.domain);
                
                % The instantiation class must recognize mtimes as a
                % functional composition.
                C.stack = @(z) A.stack(z) * B.stack(z);
                C.func = A.func * B.func;
                C.coeff = A.coeff * B.coeff;

                C.diffOrder = A.diffOrder + B.diffOrder;
            end
        end
        
        function C = plus(A, B)
            % C = A + B
            if ( isnumeric(A) )
                A = A*linBlock.eye(B.domain);
            elseif ( isnumeric(B) )
                B = B*linBlock.eye(A.domain);
            end
            dom = union(A.domain, B.domain);
            C = operatorBlock(dom);
            C.stack = @(z) A.stack(z) + B.stack(z);
            C.func = A.func + B.func;
            C.coeff = A.coeff + B.coeff;
            C.diffOrder = max(A.diffOrder, B.diffOrder);
        end
        
        function B = mpower(A, pow)
            if ( pow ~= round(pow) || pow < 0 )
                error('Power must be a positive integer.')
            end
            B = linBlock.eye(A.domain);
            for i = 1:pow
                B = B*A;
            end
        end
        
        function out = iszero(A)
            % TODO: Implement this.
            out = false;
        end
        


    end
    
end
