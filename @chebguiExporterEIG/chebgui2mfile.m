function chebgui2mfile(exporter, guifile, fid, expInfo)
%EXPORTBVP2MFILE    Export a EIG problem from CHEBGUI to a .m file.
%
%   See also: chebgui/export.
% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.

dom = expInfo.dom;
deInput = expInfo.deInput;
indVarName = expInfo.indVarName;
bcInput = expInfo.bcInput;
sigma = expInfo.sigma;
lname = expInfo.lname;
generalized = expInfo.generalized;
lhsString = expInfo.lhsString;
bcString = expInfo.bcString;
allVarNames = expInfo.allVarNames;
allVarString = expInfo.allVarString;
periodic = expInfo.periodic;
K = expInfo.K;

% Print to EIG problem
fprintf(fid, '%% Solving\n');
if ( (numel(deInput) == 1) && ~any(deInput{1} == '=') )
            fprintf(fid, '%%   %s = %s*%s\n', deInput{1}, lname, allVarString);
else
    for k = 1:numel(deInput)
        fprintf(fid, '%%   %s,\n', deInput{k});
    end
end
fprintf(fid, '%% for %s in %s', indVarName{1}, dom);
if ( ~isempty(bcInput{1}) )
    fprintf(fid,',  subject to\n%%');
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
    bcString = exporter.prettyPrintFevalString(bcString, allVarNames);
    fprintf(fid, 'N.bc = %s;\n', bcString);
end
if ( periodic )
    fprintf(fid, 'N.bc = ''periodic'';\n');
end

fprintf(fid, '\n%% Number of eigenvalue and eigenmodes to compute.\n');
fprintf(fid, 'k = %s;\n', K);

fprintf(fid, '\n%% Number of eigenvalue and eigenmodes to compute.\n');
fprintf(fid, 'k = %s;\n', K);

fprintf(fid, '\n%%%% Solve the eigenvalue problem.\n');
if ( ~generalized )
    if ( ~isempty(sigma) )
        fprintf(fid, '[V, D] = eigs(N, k, %s);\n', sigma);
    else
        fprintf(fid, '[V, D] = eigs(N, k);\n');
    end
else
    if ( ~isempty(sigma) )
        fprintf(fid, '[V, D] = eigs(N, B, k, %s);\n', sigma);
    else
        fprintf(fid, '[V, D] = eigs(N, B, k);\n');
    end
end
fprintf(fid, '\n%%%% Plot the eigenvalues.\n');
fprintf(fid, 'D = diag(D);\n');
fprintf(fid, 'figure\n');
fprintf(fid, 'plot(real(D), imag(D), ''.'', ''markersize'', 25)\n');
fprintf(fid, 'title(''Eigenvalues''); xlabel(''real''); ylabel(''imag'');\n');

if ( ischar(allVarNames) || (numel(allVarNames) == 1) )
    fprintf(fid, '\n%% Plot the eigenmodes.\n');
    fprintf(fid, 'figure\n');
    fprintf(fid, 'plot(real(V), ''linewidth'', 2);\n');
    fprintf(fid, 'title(''Eigenmodes''); xlabel(''%s''); ylabel(''%s'');', ...
        indVarName{1}, allVarString);
end

end