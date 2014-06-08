function printSetup(fid, expInfo, guifile)
%PRINTSETUP     Print commands for setting up problems
%
% Calling sequence:
%   PRINTPOSTSOLVER(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract info from the EXPINFO struct:
dom = expInfo.dom;
bcInput = expInfo.bcInput;
lname = expInfo.lname;
generalized = expInfo.generalized;
lhsString = expInfo.lhsString;
rhsString = expInfo.rhsString;
bcString = expInfo.bcString;
allVarNames = expInfo.allVarNames;
periodic = expInfo.periodic;

% Print commands for setting up problem:
fprintf(fid, '\n%%%% Define the domain we''re working on.\n');
fprintf(fid, 'dom = %s;\n', dom);
if ( ~generalized )
    fprintf(fid, ['\n%% Assign the equation to a chebop N such that' ...
        ' N(u) = %s*u.\n'], lname);
    fprintf(fid, 'N = chebop(%s, dom);\n', lhsString);
else
    fprintf(fid, ['\n%% Assign the equation to two chebops N and B such' ...
        'that N(u) = %s*B(u).\n'], lname);
    fprintf(fid, 'N = chebop(%s, dom);\n', lhsString);
    fprintf(fid, 'B = chebop(%s, dom);\n', rhsString);
end

% Make assignments for BCs.
fprintf(fid, '\n%% Assign boundary conditions to the chebop.\n');
if ( ~isempty(bcInput{1}) )
    bcString = chebguiExporter.prettyPrintFevalString(bcString, allVarNames);
    fprintf(fid, 'N.bc = %s;\n', bcString);
end
if ( periodic )
    fprintf(fid, 'N.bc = ''periodic'';\n');
end

end