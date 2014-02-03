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
    
    methods (Static = true)
        
        function F = dot(f)   % synonym for inner()
            % FUNCTIONALBLOCK.DOT Synonym for FUNCTIONALBLOCK.INNER.
            %
            % See also FUNCTIONALBLOCK.INNER.
            F = inner(f);
        end
        
        function E = eval(varargin)
%   EVAL(DOMAIN) returns a function E. The output of E(X) is an evaluator at X
%   in the domain. See also FEVAL. 
            E = @(x) functionalBlock.feval(x, varargin{:});
        end
        
        function J = jump(location,domain,order)
%JUMP      Jump at a point.
%   JUMP(LOC,DOMAIN,ORDER) returns a functional evaluating the jump (difference
%   of limit from right minus limit from left) in derivative ORDER at the point
%   LOC. 
            Er = functionalBlock.feval(location,domain,1);
            El = functionalBlock.feval(location,domain,-1);
            J = (Er-El)*operatorBlock.diff(domain,order);
        end
        
        function J = jumpAt(domain)
%JUMPAT   Jump generator. 
%   JUMP(DOMAIN) returns a callable function JMP. JMP(LOC,ORDER) calls
%   functionalBlock.jump(LOC,DOMAIN,ORDER). 
            J = @(loc,order) functionalBlock.jump(loc,domain,order);
        end
        
        function E = feval(location, varargin)
            % Use inputParser to parse the arguments to the method.
            p = inputParser;
            addRequired(p, 'location');
            pref = cheboppref;
            addOptional(p, 'domain', pref.domain, @isnumeric);
            valid = @(x) ischar(x) || isnumeric(x);
            addOptional(p, 'direction', 0, valid);
            parse(p, location, varargin{:});
            location = p.Results.location;
            domain = p.Results.domain;
            direction = p.Results.direction;
 
            % Sanity check. 
            if ( location < domain(1) ) || ( location > domain(end) )
                error('Evaluation location is not in the domain.')
            end
            
            % Convert direction argument into a number.
            if ischar(direction)
                if any( strncmpi(direction, {'left', '-'}, 1) )
                    direction = -1;
                elseif any( strncmpi(direction, {'right', '+'}, 1) )
                    direction = +1;
                else
                    error(['Direction must be ''left'', ''right'', ' ...
                        '''+'', or ''-''.'])
                end
            end
            
            % Create the FUNCTIONALBLOCK with information now available.
            E = functionalBlock(domain);
            E.stack = @(z) feval(z, location, direction);
        end
        
        function F = inner(f)
            F = functionalBlock(f.domain);
            F.stack = @(z) inner(z, f);
            F.diffOrder = 0;
        end
        
        function S = sum(domain)
            % FUNCTIONALBLOCK.SUM  Definite integration functional.
            %
            % S = FUNCTIONALBLOCK.SUM(DOMAIN) returns the definite integration
            % functional on the domain DOMAIN (i.e. the functional that maps a
            % function to its definite integral).
            pref = cheboppref;
            if ( nargin==0 )
                domain = pref.domain;
            end
            
            % Create the FUNCTIONALBLOCK with information now available.
            S = functionalBlock(domain);
            S.stack = @(z) sum(z);
            S.diffOrder = -1;
        end
        
        function Z = zero(domain)
            % FUNCTIONALBLOCK.ZERO  Zero functional.
            %
            % Z = FUNCTIONALBLOCK.ZERO(DOMAIN) returns the zero functional for
            % functions on the domain DOMAIN (i.e., the functional that maps all
            % functions to 0 on the same domain
            pref = cheboppref;
            if ( nargin == 0 )
                domain = pref.domain;
            end
            
            % Create the FUNCTIONALBLOCK with information now available.
            Z = functionalBlock(domain);
            Z.stack = @(z) zero(z);
            Z.diffOrder = 0;
        end
    end
    
end
