function bc = createBC(bcArg, ends)
%CREATEBC   Converts boundary condition syntax to CHEBFUNs so CONSTRUCTBC.m
% can deal with them.
% 
%   BC = CREATEBC(BCARG, ENDS) usually constructs a CHEBFUN for the
%   homogeneous part of the boundary conditions. If the linear constraint
%   is sufficiently compliciated then this command gives up and passes
%   resposibility to the solver. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Developer note %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This command attempts to convert the "boundary" data into CHEBFUNs. It is
% called when the boundary conditions are first assigned to the CHEBOP2 
% object. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( isa(bcArg, 'function_handle') )
    if ( nargin(bcArg) <= 1 )
        % Convert to a CHEBFUN.
        bc = chebfun(bcArg , ends);
    elseif ( nargin(bcArg) == 2 )
        % This allows for N.lbc = @(x,u) diff(u) - x + 1. 
        % We cannot convert this to a CHEBFUN so pass it on to the SOLVER.
        bc = bcArg;
    end
elseif ( isa(bcArg, 'double') )
    % This allows for N.lbc = DOUBLE. 
    bc = chebfun(bcArg, ends);
elseif ( isa(bcArg, 'chebfun') )
    % This allows for N.lbc = CHEBFUN. 
    bc = bcArg;
elseif ( isa(bcArg, 'char') )
    if ( strcmpi(bcArg, 'periodic') )
        % Pass on to SOLVER.
        bc = 'periodic'; 
    else
        error('CHEBFUN:CHEBOP2:createBC:word', ...
            'Unrecognised boundary condition string.');
    end
else
    error('CHEBFUN:CHEBOP2:createBC:type', ...
        'Unrecognised boundary condition syntax.');
end

end
