classdef linopConstraint
%LINOPCONSTRAINT Constraint class for linops.
%   This class is not intended for use by the end user. 
    
    properties
        functional    % applied to the variable to get values
        values        % constraint on the result of the functional
    end
    
    methods
        function C = linopConstraint(op, vals)
            if ( nargin == 0 )
                return
            end
            C.functional = op;
            C.values = vals;
        end
        
        function n = length(C)
            n = size(C.functional, 1);
        end
        
        function e = isempty(C)
            e = isempty(C.functional);
        end
        
        function C = append(C, op, value)
            n = length(C);
            validateattributes(op, {'linBlock', 'chebmatrix'}, {})
            if ( nargin < 3 )
                value = 0;
            end
            validateattributes(value, {'double'}, {'numel', 1})
            C.functional = [ C.functional ; op ];
            C.values(n+1, 1) = value; 
        end
        
        function C = set(C, k, op, value)
            if ( nargin < 4 )
                value = 0;
            end
            n = length(C);
            if ( k > n+1 )
                error('Requested position would leave empty constraints.')
            elseif ( k == n + 1 )
                C = append(C, op, value);
            else
                validateattributes(op, {'linBlock', 'chebmatrix'})
                validateattributes(value, {'double'}, {'numel', 1})
                C.functional{k, :} = op;
                C.values(k, 1) = value;
            end
        end

        
    end
    
end