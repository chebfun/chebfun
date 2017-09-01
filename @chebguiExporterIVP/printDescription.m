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
deInput = expInfo.deInput;
indVarNameSpace = expInfo.indVarNameSpace;
dom = expInfo.dom;
bcInput = expInfo.bcInput;

% Print a description of the problem.

% Begin by printing the differential equation:
fprintf(fid, '%% Solving\n');
for k = 1:numel(deInput)
    fprintf(fid, '%%   %s,\n', deInput{k});
end

% Print the interval:
fprintf(fid, '%% for %s in %s', indVarNameSpace, dom);

% Print the boundary conditions:
if ( ~isempty(bcInput{1}) )
    fprintf(fid, ', subject to\n%%');
    for k = 1:numel(bcInput)
        fprintf(fid, '   %s', bcInput{k});
        if ( (k ~= numel(bcInput)) && (numel(bcInput) > 1) )
            fprintf(fid, ',\n%%');
        end
    end
    fprintf(fid, '.\n');
else
    fprintf(fid, '.\n');
end

end
