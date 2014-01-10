classdef functionalBlock < linBlock
%FUNCTIONALBLOCK  Linear map of function to scalar.
%   This class is not intended to be called directly by the end user.
%
%   See also LINOP, CHEBOP, CHEBOPPREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% One of the two types of linBlock. Functionals can be composed with operators,
% added, and applied to chebfuns. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

        function C = uminus(A)
            C = functionalBlock(A.domain);
            C.stack = @(z) -A.stack(z);
        end
        
        function C = mtimes(A, B)
            % A*B 
            % If B is an operator block, the result is the composition.
            %
            % A*f 
            % If f is a chebfun, return the chebfun resulting from
            % application of A to f.
            %
            
            % Allow functional * scalar, but put the scalar first. 
            if ( isnumeric(B) && (length(B) == 1) )
                temp = A;
                A = B;
                B = temp;
            end
            
            if ( isa(B, 'chebfun') )
                % Apply functional to chebfun.
                C = toFunction(A);
                C = C(B);
            elseif ( isnumeric(A) )
                % Mulitply functional by scalar. 
                C = functionalBlock(B.domain);
                C.stack = @(z) A*B.stack(z);
                C.diffOrder = B.diffOrder;
            elseif ( isa(B, 'operatorBlock') )
                % Compose functional with operator. 
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
