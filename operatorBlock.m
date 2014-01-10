classdef operatorBlock < linBlock
%OPERATORBLOCK  Linear map of function to function.
%   This class is not intended to be called directly by the end user.
%
%   See also LINOP, CHEBOP, CHEBOPPREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% One of the two types of linBlock. Operators can be composed using *, added,
% exponentiated, and applied to chebfuns.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
    end
    
    methods
        function A = operatorBlock(domain)
            A = A@linBlock(domain);
        end
        
        function varargout = size(A, dim)
            % The dimenions of an operatorBlock are Inf-by-Inf. The syntax for
            % SIZE is the same as for a matrix. 
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
 
            % Which case?
            if ( isa(B, 'chebfun') )
                C = toFunction(A);
                C = C(B);
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
            dom = chebfun.mergeDomains(A.domain, B.domain);
            C = operatorBlock(dom);
            C.stack = @(z) A.stack(z) + B.stack(z);
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
