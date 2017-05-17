function printPostSolver(fid, expInfo)
%PRINTPOSTSOLVER   Print commands after solution has been found
%
% Calling sequence:
%   PRINTPOSTSOLVER(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract info from the EXPINFO struct:
sol = expInfo.sol;
indVarName = expInfo.indVarName;
deInput = expInfo.deInput;
s = expInfo.s;

% Plotting. Deal with scalar and system problems separately.
if ( numel(deInput) == 1 )
    % Scalar problems.
    fprintf(fid, '\n%%%% Plot the solution.\n');
    fprintf(fid, 'waterfall(%s, %s)\n', sol, ...
        indVarName{2});
    fprintf(fid, 'xlabel(''%s''), ylabel(''%s'')', ...
        indVarName{1},indVarName{2});
    
else
    % Coupled systems.
    fprintf(fid, '\n%%%% Plot the solution components.');
    for k = 1:numel(deInput)
        fprintf(fid, '\nfigure\n');
        fprintf(fid, 'waterfall(%s, %s)\n', s{k}, indVarName{2});
        fprintf(fid, 'xlabel(''%s''), ylabel(''%s''), title(''%s'')', ...
            indVarName{1},indVarName{2},s{k});
    end
    
end

end
