function printSetup(fid, expInfo, guifile)

dom = expInfo.dom;
bcInput = expInfo.bcInput;
lname = expInfo.lname;
generalized = expInfo.generalized;
lhsString = expInfo.lhsString;
bcString = expInfo.bcString;
allVarNames = expInfo.allVarNames;
periodic = expInfo.periodic;

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