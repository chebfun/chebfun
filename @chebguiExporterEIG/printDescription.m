function printDescription(fid, expInfo)
%PRINTDESCRIPTION   Print problem description when problems are exported.
%
% Calling sequence:
%   PRINTDESCRIPTION(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract info from the expInfo struct:
dom = expInfo.dom;
deInput = expInfo.deInput;
indVarName = expInfo.indVarName;
bcInput = expInfo.bcInput;
lname = expInfo.lname;
allVarString = expInfo.allVarString;
periodic = expInfo.periodic;

% Print a description of the problem.

% Begin by printing the differential equation:
fprintf(fid, '%% Solving\n');
if ( (numel(deInput) == 1) && ~any(deInput{1} == '=') )
    fprintf(fid, '%%   %s = %s*%s\n', deInput{1}, lname, allVarString);
else
    for k = 1:numel(deInput)
        fprintf(fid, '%%   %s,\n', deInput{k});
    end
end

% Print the interval that x lives on:
fprintf(fid, '%% for %s in %s', indVarName{1}, dom);

% Print the boundary conditions:
if ( ~isempty(bcInput{1}) )
    % Non-periodic conditions:
    fprintf(fid,', subject to\n%%');
    for k = 1:numel(bcInput)
        fprintf(fid, '   %s', bcInput{k});
        if ( (k ~= numel(bcInput)) && (numel(bcInput) > 1) )
            fprintf(fid, ',\n%%');
        end
    end
    fprintf(fid, '.\n');
elseif ( periodic )
    fprintf(fid, ', subject to periodic boundary conditions.\n\n');
else
    fprintf(fid, '.\n');
end

end
