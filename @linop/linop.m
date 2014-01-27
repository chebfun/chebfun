classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) linop < chebmatrix
%LINOP     Linear operator with boundary and side conditions.
%   A linop is a chebmatrix plus side constraints, representing boundary
%   conditions, for example. 
%
%   If A is a chebmatrix, then L = LINOP(A) converts it to a linop. Use
%   LINOP.ADDBC to attach boundary and side conditions to L. 
%
%   If A,B,C,... are linops, then concatenation such as L = [ A, B, C ]
%   creates a linop with no side conditions. The sizes and domains of the
%   individual parts must be compatible. 
%
%   There is no need to directly construct a LINOP to solve a differential
%   equation. Instead, use a CHEBOP, which creates the needed linops
%   automatically.
%
%   Example:
%     d = [-2 2];                       % function domain
%     I = operatorBlock.eye(d);         % identity
%     D = operatorBlock.diff(d);        % differentiation
%     x = chebfun(@(x) x, d);           % the variable x on d
%     M = operatorBlock.mult(x.^2);     % multiplication operator
%     S = functionalBlock.sum(d);       % integration functional 
%     E = functionalBlock.eval(d);      % evaluation functional generator
% 
%     u = [ exp(x); pi; sin(x) ];       %  function; scalar; function
%     A = [ I+D, abs(x), M;
%           S, 0, E(2);  
%           I, x.^2, D ];
% 
%     A = linop(A);
%     z = functionalBlock.zeros(d);
%     A = addbc( A, [E(-2), 0, z], 0 );         % set u{1}(-2) = 0
%     A = addbc( A, [E(2)*D, 0, E(2)], 0 );     % set u{1}'(2) + u{3}(2) = 0
%
%     f = [ x; 1; 0*x ];
%     u = A\f;
%     plot( chebfun(u) )
%
%   See also CHEBOPPREF, CHEBOP, CHEBMATRIX, LINOP.ADDBC.
    
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% The LINOP class is essentially a CHEBMATRIX with two additional properties:
% linop.constraint and linop.continuity which allow imposing conditions on the
% solution (e.g. boundary conditions) and continuity conditions. Both
% linop.constraint and linop.continuity are of the type LINOPCONSTRAINT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        constraint = linopConstraint()
        continuity = linopConstraint()
    end
    
    methods
        function L = linop(M)
            L = L@chebmatrix(M);
        end
    end
    
    methods ( Static )
        function [Z, I, D, C] = primitiveOperators(domain)
        %LINOP.PRIMITIVEOPERATORS Frequently used operator blocks.
        %   [Z, I, D, C] = LINOP.PRIMITIVEOPERATORS(DOMAIN) returns a zero,
        %   identity, differentiation, and indefinite integration operator,
        %   respecitvely, on the given domain.
        %
        %   Note that the 'linop' is just a convenient name. The outputs are
        %   *not* linops.
        %
        %   See also OPERATORBLOCK.
            
            Z = operatorBlock.zeros(domain);
            I = operatorBlock.eye(domain);
            D = operatorBlock.diff(domain,1);
            C = operatorBlock.cumsum(domain,1);
        end

        function [z, E, s] = primitiveFunctionals(domain)
        %LINOP.PRIMITIVEFUNCTIONALS Frequently used functional blocks.
        %   [Z, E, S] = LINOP.PRIMITIVEFUNCTIONALS(DOMAIN) returns frequently
        %   used functional blocks on the given domain.
        %
        %   Z is the zero functional. E(X) is an evaluation functional at the
        %   point X. E(X, DIR) adds a direction (+1/-1) to the evaluation. S is
        %   a definite integration functional.
        %
        %   Note that the 'linop' is just a convenient name. The outputs are
        %   *not* linops.
        %
        %   See also FUNCTIONALBLOCK.
            
            z = functionalBlock.zero(domain);
            E = functionalBlock.eval(domain);
            s = functionalBlock.sum(domain);
        end
    end

end
