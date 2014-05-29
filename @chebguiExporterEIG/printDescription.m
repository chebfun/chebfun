function printDescription(fid, expInfo)
%PRINTDESCRIPTION   Print problem description when problems are exported.

% Extract info from the expInfo struct:
deInput = expInfo.deInput;
indVarNameSpace = expInfo.indVarNameSpace;
dom = expInfo.dom;
bcInput = expInfo.bcInput;
periodic = expInfo.periodic;

%% Print a description of the BVP:

% Print the differential equation:
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
elseif periodic
    fprintf(fid, ', subject to periodic boundary conditions.\n\n');
else
    fprintf(fid, '.\n');
end

end