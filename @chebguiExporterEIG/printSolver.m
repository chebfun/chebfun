function printSolver(fid, expInfo)
sigma = expInfo.sigma;
generalized = expInfo.generalized;

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
end