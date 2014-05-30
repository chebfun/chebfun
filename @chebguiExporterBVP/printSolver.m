function printSolver(fid, expInfo)
%PRINTSOLVER    Print the solution step when exporting
fprintf(fid,'\n%%%% Solve the problem!');
fprintf(fid, ['\n%% Here, we call the solvebvp() method ' ...
    '(which offers the same functionality \n%% as nonlinear '...
    'backslash, but with more customizability).\n']);
fprintf(fid, 'u = solvebvp(N, rhs, options);\n');
end