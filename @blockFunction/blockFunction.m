classdef blockFunction
%BLOCKFUNCTION   Convert linear operator to callable function.
%   This class is not intended to be called directly by the end user.
%
%   See also LINOP, CHEBOP, CHEBOPPREF.
    
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% This class converts a linBlock object into a callable function suitable
% for application to a chebfun.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties ( Access=public )
        % This property is assigned the callable function that does the
        % correct operation when called on CHEBFUN objects.
        func = [];
        domain;
    end
    
    methods
        
        function A = blockFunction(varargin)
            % When called with no arguments, the returned object causes the
            % block's stack to be evaluated with these methods to produce
            % coefficients.
            if ( isempty(varargin{1}) )
                pref = cheboppref;
                A.domain = pref.domain;  % not clear this is ever used...
                return
                
                
           % Calling the constructor with a linBlock argument initiates the
           % process of evaluating the stack with a dummy object of this class.
           elseif ( isa(varargin{1}, 'linBlock') )
                % Convert the given linBlock to its function form by
                % evaluating its stack.
                L = varargin{1};
                dummy = blockFunction([]);
                dummy.domain = L.domain;
                A = L.stack( dummy );
                
                
            % If the constructor is called with data, just make a regular object
            % out of it. 
            else
                 A.func = varargin{1};
            end
        end
        
        
        function C = mtimes(A, B)
            % Interpret as composition.
            C = blockFunction( @(z) A.func(B.func(z)) );
        end
        
        function C = plus(A, B)
            C = blockFunction( @(z) A.func(z) + B.func(z) );
        end
        
        function C = uminus(A)
            C = blockFunction( @(z) -A.func(z) );
        end
        
        function A = uplus(A)
        end
        
    end
    
    methods
        function D = diff(A, m)
            D = blockFunction( @(z) diff(z,m) );
        end
        
        function C = cumsum(A, m)
            C = blockFunction( @(u) cumsum(u,m) );
        end
        
        function I = eye(A)
            I = blockFunction( @(z) z );
        end
        
        function Z = zeros(A)
            Z = blockFunction( @(z) chebfun(0,A.domain) );
        end
        
        function z = zero(A)
            z = blockFunction( @(u) 0 );
        end
        
        function F = mult(A, f)
            F = blockFunction( @(z) times(f,z) );
        end
        
        function S = sum(A)
            S = blockFunction( @(z) sum(z) );
        end
        
        function E = feval(A, location, direction)
            if ( direction < 0 )
                E = blockFunction( @(u) feval(u, location, -1) );
            elseif ( direction > 0 )
                E = blockFunction( @(u) feval(u, location, 1) );
            else
                E = blockFunction( @(u) feval(u, location) );
            end
        end
        
        function F = inner(A, f)
            F = blockFunction( @(z) mtimes(f', z) );
        end
        
        
    end
    
end
