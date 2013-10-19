classdef functionalBlock < linBlock
    
    properties
    end
    
    methods
        function A = functionalBlock(domain)
            A = A@linBlock(domain);
        end
        
        function varargout = size(A, dim)
            % S = SIZE(A)
            % [M, N] = SIZE(A)
            % P = SIZE(A, K)

            m = [1, Inf];
            if nargin > 1
                varargout = {m(dim)};
            elseif nargout <= 1
                varargout = {m};
            else
                varargout = {m(1) m(2)};
            end
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
            
            if ( isnumeric(B) && length(B) == 1 )
                C = mtimes(B, A);
            end
            
            if ( isa(B, 'chebfun') )
                C = op(A);
                C = C(B);
            elseif ( isnumeric(A) )
                C = functionalBlock(A.domain);
                C.stack = A*B.stack(z);
                C.diffOrder = B.diffOrder;
            elseif ( isa(B, 'operatorBlock') )
                C = functionalBlock(A.domain);
                C.stack = @(z) A.stack(z) * B.stack(z);
                C.diffOrder = A.diffOrder + B.diffOrder;
            else 
                error('Unrecognized operand types.')
            end
        end
        
        function C = plus(A, B)
            % C = A + B
            C = functionalBlock(A.domain);
            C.stack = @(z) A.stack(z) + B.stack(z);
            C.diffOrder = max(A.diffOrder, B.diffOrder);
        end        
        

    end
    
end
