function display( N )
%DISPLAY   Print to command line for CHEBOP2.
% DISPLAY is called automatically when a statement that results in a
% CHEBOP2 output is not terminated with a semicolon.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

loose = ~isequal(get(0,'FormatSpacing'),'compact');

if ( loose ) 
    fprintf('\n') 
end

% compact version
if ( isempty( N ) )
    fprintf('    empty chebop2\n\n')
    return
end

% Get information that we want to display
dom = N.domain;                           % Domain
op = N.op; 

% Display the information: 

disp('chebop2 object: (partial differential equation)')
fprintf('         domain               operator  \n');

fprintf(['[%4.2g,%4.2g] x [%4.2g,%4.2g]    ' func2str(op)  '\n'], dom(1),...
                                          dom(2), dom(3), dom(4));
     
PrintBoundaryConditions(N.lbc, 'left');  
PrintBoundaryConditions(N.rbc, 'right');
PrintBoundaryConditions(N.ubc, 'top');
PrintBoundaryConditions(N.dbc, 'bottom');

if ( loose ) 
    fprintf('\n');
end

end

function PrintBoundaryConditions(BC, string)
% Display the boundary conditions: 

    exists = ~isempty(BC);
    if ( exists )
        if ( isa(BC, 'double') )
            fprintf(['boundary ' string ' = %u'], BC);
        elseif ( isa(BC, 'chebfun') )
            % Try to convert to a double: 
            n = length(BC); 
            if ( n == 1 )
                val = BC.values{:};
                PrintBoundaryConditions(val, string)
                return
            end
            fprintf(['boundary ' string ' = function']);
        elseif ( isa(BC, 'function_handle') )
            fprintf(['boundary ' string ' = ' func2str(BC)]);
        end
        fprintf('\n');
    end
end