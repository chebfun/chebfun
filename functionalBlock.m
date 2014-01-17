classdef functionalBlock < linBlock
%FUNCTIONALBLOCK  Linear map of function to scalar.
%   This class is not intended to be called directly by the end user.
%
%   See also LINBLOCK, LINOP, CHEBOP, CHEBOPPREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% One of the two concrete implementations of the abstract type LINBLOCK.
% Functionals can be composed with operators, added, and applied to CHEBFUN
% objects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    methods
        function A = functionalBlock(domain)
            % FUNCTIONALBLOCK constructor, simply calls the LINBLOCK
            % constructor.
            A = A@linBlock(domain);
        end
        
        function varargout = size(A, dim) %#ok<INUSL>
            % SIZE  Size of a FUNCTIONALBLOCK.
            % The commands
            %   S = SIZE(A)
            %   [M, N] = SIZE(A)
            %   P = SIZE(A, K)
            % return the results expected from standard MATLAB syntax.
            
            % A FUNCTIONALBLOCK is always of dimensions 1 x Inf.
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
            % Unary minus of a FUNCTIONALBLOCK
            C = functionalBlock(A.domain);
            C.stack = @(z) -A.stack(z);
        end
        
        function C = mtimes(A, B)
            % *    Functional composition, multiplication or application.
            %
            % C = A*B, where A is a FUNCTIONALBLOCK and B is an OPERATORBLOCK,
            % returns the FUNCTIONALBLOCK C which is the the result of composing
            % the operators A and B.
            %
            % C = A*B, or C = B*A,  where A is a FUNCTIONALBLOCK and B is a
            % scalar, returns the FUNCTIONALBLOCK C which is the the result of
            % multiplying A with B.
            %
            % C = A*F, where A is a FUNCTIONALBLOCK and F is a CHEBFUN, returns
            % the CHEBFUN C which is the the result of applying A to F.

            
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
            % C = A + B     Addition of two FUNCTIONALBLOCK objects A and B.
            C = functionalBlock(A.domain);
            C.stack = @(z) A.stack(z) + B.stack(z);
            C.diffOrder = max(A.diffOrder, B.diffOrder);
        end        
        

    end
    
end
