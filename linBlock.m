classdef linBlock
%LINBLOCK   Linear operator on a single function.
%   This class is not intended to be called directly by the end user.
%
%   See also LINOP, CHEBOP, CHEBOPPREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% A LINBLOCK is an abstract representation of a linear operator on a single
% function defined on a fixed domain. Its main purpose is to maintain an
% unevaluated stack of algebraic steps operating on predefined building blocks.
% It also keeps track of the differential order of the operator.
%
% The reason for the 'Block' name is that these objects can be a block entry in
% a chebmatrix, for operators that apply to a mixture of multiple functions and
% scalars. The preferred usage externally is to create a 1x1 chebmatrix or linop
% rather than to manipulate the blocks themselves.
%
% The stack is used three different ways within the system:
%       (1) Replace the  building blocks with appropriate matrices, and then
%           apply the algebra. (I.e., collocation.) 
%       (2) Derive the coefficients of each power of the derivative. 
%           (See the method toCoeff.) 
%       (3) Convert to a callable function that can be applied directly to a
%           CHEBFUN. (See the method toFunction.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties
        % The domain of CHEBFUN objects that are operated upon.
        domain = [];
        
        % The delayed evaluation stack. Its argument is an empty
        % instance of the class that determines how the operator is
        % discretized.
        stack = [];
                
        % Used to track differential order.
        diffOrder = 0;        
    end
    
    methods
        function A = linBlock(varargin)
            % A = LINBLOCK()        Use the default preference domain. 
            % A = LINBLOCK(B)       Self-return for the same type.
            % A = LINBLOCK(DOMAIN)  Null object on the domain.
            
            if ( nargin == 0 )
                p = chebpopref;
                A.domain = p.domain;
            elseif ( nargin == 1 && isa(varargin{1}, 'linBlock') )
                A = varargin{1};
                return
            else
                A.domain = varargin{1};
            end
        end
        
        function f = toFunction(A)
            % Convert the LINBLOCK to a callable anonymous function that can be
            % applied to a CHEBFUN.
            B = blockFunction(A);
            f = B.func;
        end
        
        function c = toCoeff(A)
            % Convert the LINBLOCK to a CHEBFUN of coefficients multiplying the
            % powers of the derivative.
            B = blockCoeff(A);
            c = chebmatrix( B.coeffs );
        end
        
        function C = minus(A, B)
            C = A + (-B);
        end
               
        function C = horzcat(varargin)
            kill = cellfun(@isempty, varargin);   % skip empties         
            C = chebmatrix( varargin(~kill) );
        end
        
        function C = vertcat(varargin)
            kill = cellfun(@isempty, varargin);   % skip empties         
            C = chebmatrix( varargin(~kill)' );
        end
                                        
    end
        
    methods (Static)
        function D = diff(varargin)
            % LINBLOCK.DIFF  Differentiation operator.
            %
            % D = LINBLOCK.DIFF returns the first-order differentation operator
            % D for functions defined on [-1, 1].
            %
            % D = LINBLOCK.DIFF(DOMAIN) returns the first-order differentation
            % operator D which applies to functions defined on DOMAIN, which may
            % include breakpoints.
            %
            % D = LINBLOCK.DIFF(DOMAIN, M) is the mth order derivative.
            
            % Use inputParser to parse the arguments to the method.
            p = inputParser;
            pref = cheboppref;
            addOptional(p, 'domain', pref.domain, @isnumeric);
            mcheck = @(m) validateattributes(m, ...
                {'numeric'}, {'scalar', 'nonnegative', 'integer'});
            addOptional(p, 'm', 1, mcheck);
            parse(p, varargin{:});
            dom = p.Results.domain;
            m = p.Results.m;
            
            % Create the LINBLOCK with information now available. This LINBLOCK
            % will be off the concrete type OPERATORBLOCK.
            D = operatorBlock(dom);
            D.stack = @(z) diff(z, m);
            D.diffOrder = m;
        end
        
        function C = cumsum(varargin)
            % LINBLOCK.CUMSUM  Antiderivative operator.
            %
            % C = LINBLOCK.CUMSUM returns the first-order antiderivative operator C
            % for functions defined on [-1, 1]. The result of applying the
            % operator is defined uniquely by having value zero at the left
            % endpoint.
            %
            % C = LINBLOCK.CUMSUM(DOMAIN) returns the first-order antiderivative
            % operator C which applies to functions defined on DOMAIN, which may
            % include breakpoints.
            %
            % C = LINBLOCK.CUMSUM(DOMAIN, M) is the mth-repeated antiderivative.
            
            % Use inputParser to parse the arguments to the method.
            p = inputParser;
            pref = cheboppref;
            addOptional(p, 'domain', pref.domain, @isnumeric);
            mcheck = @(m) validateattributes(m, {'numeric'}, ...
                {'scalar', 'nonnegative', 'integer'});
            addOptional(p, 'm', 1, mcheck);
            parse(p, varargin{:});
            dom = p.Results.domain;
            m = p.Results.m;
            
            % Create the LINBLOCK with information now available. This LINBLOCK
            % will be off the concrete type OPERATORBLOCK.
            C = operatorBlock(dom);
            C.stack = @(z) cumsum(z, m);
            C.diffOrder = -m;
        end
        
        function I = eye(domain)
            % LINBLOCK.EYE  Identity operator.
            %
            % I = LINBLOCK.EYE(DOMAIN) returns the identity operator for
            % functions on the domain DOMAIN.
            pref = cheboppref;
            if ( nargin == 0 )
                domain = pref.domain; 
            end
            
            % Create the LINBLOCK with information now available. This LINBLOCK
            % will be off the concrete type OPERATORBLOCK.
            I = operatorBlock(domain);
            I.stack = @(z) eye(z);
            I.diffOrder = 0;
        end
        
        function Z = zeros(domain)
            % LINBLOCK.ZEROS  Zero operator.
            %
            % Z = LINBLOCK.ZEROS(DOMAIN) returns the zero operator for functions
            % on the domain DOMAIN (i.e., the operator that maps all functions
            % to the zero function on the DOMAIN).
            pref = cheboppref;
            if ( nargin == 0 )
                domain = pref.domain;
            end
            
            % Create the LINBLOCK with information now available. This LINBLOCK
            % will be off the concrete type OPERATORBLOCK.
            Z = operatorBlock(domain);
            Z.stack = @(z) zeros(z);
            Z.diffOrder = 0;
        end
        
        function Z = zero(domain)
            % LINBLOCK.ZERO  Zero functional.
            %
            % Z = LINBLOCK.ZERO(DOMAIN) returns the zero functional for
            % functions on the domain DOMAIN (i.e., the functional that maps all
            % functions to 0on the same domain
            pref = cheboppref;
            if ( nargin == 0 )
                domain = pref.domain;
            end
            
            % Create the LINBLOCK with information now available. This LINBLOCK
            % will be off the concrete type FUNCTIONALBLOCK.
            Z = functionalBlock(domain);
            Z.stack = @(z) zero(z);
            Z.diffOrder = 0;
        end
        
        function M = mult(u)
            % LINBLOCK.MULT  Multiplication operator.
            %
            % M = LINBLOCK.MULT(U) returns the multiplication operator from the
            % CHEBFUN U, i.e. the operator that maps a CHEBFUN f(x) to u(x)f(x).
            
            % Create the LINBLOCK with information now available. This LINBLOCK
            % will be off the concrete type OPERATORBLOCK.
            M = operatorBlock(u.domain);
            M.stack = @(z) mult(z, u);
            M.diffOrder = 0;
        end
        
        function S = sum(domain)
            % LINBLOCK.SUM  Definite integration functional.
            %
            % S = LINBLOCK.SUM(DOMAIN) returns the definite integration
            % functional on the domain DOMAIN (i.e. the functional that maps a
            % function to its definite integral).
            pref = cheboppref;
            if ( nargin==0 )
                domain = pref.domain;
            end
                        
            % Create the LINBLOCK with information now available. This LINBLOCK
            % will be off the concrete type FUNCTIONALBLOCK.
            S = functionalBlock(domain);
            S.stack = @(z) sum(z);
            S.diffOrder = -1;
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
            
            % Create the LINBLOCK with information now available. This LINBLOCK
            % will be off the concrete type FUNCTIONALBLOCK.
            E = functionalBlock(domain);
            E.stack = @(z) feval(z, location, direction);
        end
        
        function E = eval(varargin)
            % EVAL(DOMAIN) returns a function E. The output of E(X) is an
            % evaluator at X in the domain.
            E = @(x) linBlock.feval(x, varargin{:});
        end
        
        function F = inner(f)
            F = functionalBlock(f.domain);
            F.stack = @(z) inner(z, f);
            F.diffOrder = 0;
        end

        function F = dot(f)   % synonym for inner()
            % LINBLOCK.DOT Synonym for LINBLOCK.INNER.
            %
            % See also LINBLOCK.INNER.
            F = inner(f);
        end
        
        function F = fred(domain, kernel, varargin)
            F = operatorBlock(domain);
            F.stack = @(z) fred(z, kernel, varargin{:});
            F.diffOrder = 0;
        end

        function V = volt(domain, kernel, varargin)
            V = operatorBlock(domain);
            V.stack = @(z) volt(z, kernel, varargin{:});
            V.diffOrder = 0;
        end

        
    end
    

    
end