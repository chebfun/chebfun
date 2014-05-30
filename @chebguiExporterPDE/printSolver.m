function printSolver(fid, expInfo)
indVarName = expInfo.indVarName;
sol = expInfo.sol;
sol0 = expInfo.sol0;
deInput = expInfo.deInput;
s = expInfo.s;

fprintf(fid, ['\n%% Solve the problem using pde15s.\n']);
fprintf(fid, '[%s, %s] = pde15s(pdefun, %s, %s, bc, opts);\n', indVarName{2}, ...
    sol, indVarName{2}, sol0);

% Conver sol to variable names
if ( numel(deInput) > 1 )
    fprintf(fid, '\n%% Recover variable names.\n');
    for k = 1:numel(s)
        fprintf(fid, '%s = %s{%d};\n', s{k}, sol, k);
    end
end
end