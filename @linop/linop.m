classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) linop < chebmatrix
%LINOP     Linear operator with boudnary and side conditions.
%   A linop is a chebmatrix plus side constraints, representing boundary
%   conditions, for example. 
%
%   If A is a chebmatrix, then L = LINOP(A) converts it to a linop. Use
%   ADDBC to attach boundary and side conditions to L. 
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
%     d = [-2 2];   % function domain
%     I = operatorBlock.eye(d);   % identity
%     D = operatorBlock.diff(d);  % differentiation
%     x = chebfun(@(x)x,d);   % the variable x on d
%     M = operatorBlock.mult(x.^2);   % multiplication operator
%     S = functionalBlock.sum(d);   % integration functional 
%     E = functionalBlock.eval(d);  % evaluation functional generator
% 
%     u = [ exp(x); pi; sin(x) ];   %  function; scalar; function
%     A = [ I+D, abs(x), M;
%           S, 0, E(2);  
%           I, x.^2, D ];
% 
%     A = linop(A);
%     z = functionalBlock.zero(d);
%     A = addbc( A, [E(-2),0,z], 0 );   % set u{1}(-2) = 0
%     A = addbc( A, [E(2)*D,0,E(2)], 0 );  % set u{1}'(2) + u{3}(2) = 0
%
%     f = [ x; 1; 0*x ];
%     u = A\f;
%     plot( chebfun(u) )
%
%   See also CHEBOPPREF, CHEBOP, CHEBMATRIX, LINOP.ADDBC.
    
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
         
    properties
        constraint = linopConstraint()
        continuity = linopConstraint()
    end
    
    methods
        function L = linop(M)
            if ( ~strcmp(class(M),'chebmatrix') )
                error('Input must be a chebmatrix.')
            end
            L = L@chebmatrix(M);
        end
        
    end
<<<<<<< Updated upstream
    
end
=======
        
end
>>>>>>> Stashed changes
