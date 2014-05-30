function printPostSolver(fid, expInfo)
%PRINTPOSTSOLVER    Print commands after solution has been found

% Extract information from the EXPINFO struct
allVarNames = expInfo.allVarNames;
indVarNameSpace = expInfo.indVarNameSpace;

% Print commands that will create a plot of the solution obtained:
fprintf(fid, '\n%%%% Create a plot of the solution.\n');

fprintf(fid, ['figure\nplot(u,''LineWidth'',2)\n', ...
    'title(''Final solution''), xlabel(''%s'')'], indVarNameSpace);
if ( numel(allVarNames) == 1 )
    % Scalar problem:
    fprintf(fid, ', ylabel(''%s'')', allVarNames{:});
else
    % Coupled system. Create a legend.
    leg = '';
    for k = 1:numel(allVarNames)-1
        leg = [leg '''' allVarNames{k} '''' ','];
    end
    leg = [leg '''' allVarNames{k+1} ''''];
    fprintf(fid, ', legend(%s)\n', leg);
end


end
