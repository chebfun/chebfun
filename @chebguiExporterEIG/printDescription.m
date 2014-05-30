function printDescription(fid, expInfo)
%PRINTDESCRIPTION   Print problem description when problems are exported.

% Extract info from the expInfo struct:
dom = expInfo.dom;
deInput = expInfo.deInput;
indVarName = expInfo.indVarName;
bcInput = expInfo.bcInput;
lname = expInfo.lname;
allVarString = expInfo.allVarString;
periodic = expInfo.periodic;

%% Print a description of the eigenvalue problem:

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

end