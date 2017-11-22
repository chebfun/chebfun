classdef linopConstraint
%LINOPCONSTRAINT Class for linop constraints.
%
%   A linop operates on a set of chebfuns and scalars (the "variables" of the
%   linop). An instance of this class may be assigned to the 'constraint'
%   property of a linop in order to impose a constraint on that variable (for
%   purposes of solving a linear system or eigenvalue problem).
% 
%   Each LINOPCONSTRAINT object has a 'functional' property, which is a
%   chebmatrix of vertically concatenated functionalBlocks, and a 'values'
%   property, a vector of the same length. If u is the variable of the
%   linop, then the constraint imposed by object C is:
%
%        C.functional * u = C.values.
%
%   The preferred way to add contraints to a linop is through LINOP.ADDBC,
%   not through direct use of this class.
%
% See also LINOP, LINOP.ADDBC.

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        functional    % applied to the variable to get values
        values        % constraint on the result of the functional
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        % Input a functional and a value to create a constraint, or create
        % an empty constraint if no inputs.
        function C = linopConstraint(op, vals)
            if ( nargin == 0 )
                return
            end
            C.functional = op;
            C.values = vals;
        end
        
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function n = length(C)
        %LENGTH    Number of constraints in the object.            
            n = size(C.functional, 1);
        end
        
        function e = isempty(C)
        %ISEMPTY   True if no constraints in the object.            
            e = isempty(C.functional);
        end
        
        function C = append(C, func, value)
        %APPEND    Insert an additional constraint.
        %   C = APPEND(C, FUNC, VAL) appends the constraint FUNC*u=VAL to the
        %   current list. If VAL is omitted, it defaults to zero.
        
            % Appending empty does nothing:
            if ( isempty(func) )
                return
            end                
            
            % Check if func is of an allowed class.
            validateattributes(func, {'linBlock', 'chebmatrix'}, {})
            
            % Check the value.
            if ( nargin < 3 )
                value = 0;
            end
            validateattributes(value, {'double'}, {})
            
            C.functional = [ C.functional ; func ];
            C.values = [ C.values ; value ];
        end
        
        function C = uminus(C)
        %-   Negate the VALUES property of a LINOPCONSTRAINT.
        %   This is useful at the CHEBOP level, where we need different signs
        %   for the boundary conditions of a LINOP when we call LINOP(N) where N
        %   is a CHEBOP, compared to what we want when we call LINEARIZE() from
        %   within a Newton iteration. This is because when problems are solved
        %   with LINOP backslash, the solution to the problem is the output
        %   itself, while in a Newton iteration, we have to add the output of
        %   the LINOP solution to the current guess.
            C.values = -C.values;
        end
        
    end
    
end
