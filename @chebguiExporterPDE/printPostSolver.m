function printPostSolver(fid, expInfo)
sol = expInfo.sol;
indVarName = expInfo.indVarName;
deInput = expInfo.deInput;
s = expInfo.s;

% plotting
if ( numel(deInput) == 1 )
    fprintf(fid, '\n%% Create plot of the solution.\n');
    %     fprintf(fid,'surf(%s,t,''facecolor'',''interp'')\n',sol);
    fprintf(fid, 'waterfall(%s,%s,''simple'',''linewidth'',2)\n', sol, ...
        indVarName{2});
else
    fprintf(fid, '\n%% Create plots of the solutions.\n');
    M = numel(deInput);
    for k = 1:numel(deInput)
        fprintf(fid, 'subplot(1,%d,%d)\n', M, k);
        fprintf(fid, 'waterfall(%s,%s,''linewidth'',2)\n', s{k}, indVarName{2});
        fprintf(fid, 'xlabel(''%s''), ylabel(''%s''), title(''%s'')\n', ...
            indVarName{1},indVarName{2},s{k});
    end
end

end