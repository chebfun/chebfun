function printOptions(fid, expInfo)
%PRINTOPTIONS   Print problem options when they are exported.
%
% Calling sequence:
%   PRINTOPTIONS(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract info from the expInfo struct:
tol = expInfo.tol;
dampingOn = expInfo.dampingOn;
discretization = expInfo.discretization;
timeSteppingSolver = expInfo.timeSteppingSolver;
ivpSolver = expInfo.ivpSolver;
plotting = expInfo.plotting;

% Set up preferences
fprintf(fid, '\n%%%% Setup preferences for solving the problem.\n');
fprintf(fid, '%% Create a CHEBOPPREF object for passing preferences.\n');
fprintf(fid, '%% (See ''help cheboppref'' for more possible options.)\n');
   
fprintf(fid, 'options = cheboppref();\n\n');

% Specify which IVP solver we want to use:
fprintf(fid, '%% Specify the IVP solver to use. Possible options are:\n');
fprintf(fid, '%%   Time-stepping solvers:\n');
fprintf(fid, '%%     ''ode113'' (default), ''ode15s'' or ''ode45''.\n');
fprintf(fid, '%%   Global methods:\n');
fprintf(fid, '%%     ''values'' or ''coefficients''.\n');
fprintf(fid, 'options.ivpSolver = ''%s'';\n', ivpSolver);

% We'll need different preferences depending on whether we're applying a global
% or a time-stepping solver:
if ( timeSteppingSolver )
    % Do nothing else (for now).
else
    % Always show iteration information:
    fprintf(fid, 'options.display = ''iter'';\n');
    
    % Specify tolerance:
    if ( ~isempty(tol) )
        fprintf(fid, '\n%% Option for tolerance.\n');
        fprintf(fid, 'options.bvpTol = %s;\n', tol);
    end
    
    % Option for damping:
    fprintf(fid, '\n%% Option for damping.\n');
    if ( strcmp(dampingOn, '1') )
        fprintf(fid, 'options.damping = true;\n');
    else
        fprintf(fid, 'options.damping = false;\n');
    end
    
    % Option for discretization:
    fprintf(fid, ['\n%% Option for discretization (either ''values'' ' ...
        'or ''coeffs'').\n']);
    fprintf(fid, 'options.discretization = ''%s'';\n', ivpSolver);
    
    % Plot during Newton iteration?
    if ( ~strcmp(plotting, 'off') )
        fprintf(fid, ['\n%% Option for determining how long each Newton step ' ...
            'is shown.\n']);
        fprintf(fid, 'options.plotting = %s;\n', plotting);
    else
        fprintf(fid, ['\n%% Option for determining how long each Newton step ' ...
            'is shown.\n']);
        fprintf(fid,'options.plotting = ''off'';\n');
    end
end
end
