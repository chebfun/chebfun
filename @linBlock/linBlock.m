classdef (InferiorClasses = {?chebfun}) linBlock
%LINBLOCK   Linear operator on a single function.
%   This class is not intended to be called directly by the end user.
%
% See also LINOP, CHEBOP, CHEBOPPREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % The domain of CHEBFUN objects that are operated upon.
        domain = [];
        
        % The delayed evaluation stack. Its argument is an empty
        % instance of the class that determines how the operator is
        % discretized.
        stack = [];
                
        % Used to track differential order.
        diffOrder = 0;
        
        % Is the operator the zero operator or the zero functional? Usually not,
        % so default value is set to 0: 
        iszero = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function A = linBlock(varargin)
            % A = LINBLOCK()        Use the default preference domain. 
            % A = LINBLOCK(B)       Self-return for the same type.
            % A = LINBLOCK(DOMAIN)  Null object on the domain.
            
            if ( nargin == 0 )
                p = cheboppref();
                A.domain = p.domain;
            elseif ( (nargin == 1) && isa(varargin{1}, 'linBlock') )
                A = varargin{1};
                return
            else
                A.domain = varargin{1};
            end
            
        end
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
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
        
        function spy(A, varargin)
           % The CHEBMATRIX class implements a SPY method, which is precisely
           % what we want:
           spy(chebmatrix({A}), varargin{:});
        end
        
        function C = vertcat(varargin)
            kill = cellfun(@isempty, varargin);   % skip empties         
            C = chebmatrix( varargin(~kill)' );
        end
                                        
    end
    
end
