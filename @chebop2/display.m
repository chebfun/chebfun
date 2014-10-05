function display(N)
%DISPLAY   Print to command line for CHEBOP2.
% DISPLAY is called automatically when a statement that results in a
% CHEBOP2 output is not terminated with a semicolon.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

loose = ~isequal(get(0, 'FormatSpacing'), 'compact');

if ( loose ) 
    fprintf('\n') 
end

% Compact version.
if ( isempty(N) )
    fprintf('    empty chebop2\n\n')
    return
end

% Get information that we want to display.
dom = N.domain;                           
op = N.op; 

% Display the information.
fprintf('chebop2 of rank %u: (linear bivariate operator)\n',rank(N))
fprintf('         domain               operator  \n');

fprintf(['[%3.2g,%3.2g] x [%3.2g,%3.2g]  ' func2str(op)  '\n'], dom(1),...
                                          dom(2), dom(3), dom(4));
printBoundaryConditions(N.lbc, 'left');  
printBoundaryConditions(N.rbc, 'right');
printBoundaryConditions(N.ubc, 'top');
printBoundaryConditions(N.dbc, 'bottom');

if ( loose ) 
    fprintf('\n');
end

end

function printBoundaryConditions(BC, string)
% Display the boundary conditions.

    exists = ~isempty(BC);
    if ( exists )
        if ( isa(BC, 'double') )
            fprintf([string ' boundary = %u'], BC);
        elseif ( isa(BC, 'chebfun') )
            % Try to convert to a double: 
            n = length(BC); 
            if ( n == 1 )
                val = BC.values;
                printBoundaryConditions(val, string)
                return
            end
            fprintf([string ' boundary = function']);
        elseif ( isa(BC, 'function_handle') )
            fprintf([string ' boundary = ' func2str(BC)]);
        end
        fprintf('\n');
    end
    
end