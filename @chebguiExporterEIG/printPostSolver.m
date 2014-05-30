function printPostSolver(fid, expInfo)

allVarNames = expInfo.allVarNames;
indVarName = expInfo.indVarName;
allVarString = expInfo.allVarString;

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