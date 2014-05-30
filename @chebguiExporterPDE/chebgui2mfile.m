function chebgui2mfile(exporter, guifile, fid, expInfo)
%EXPORTBVP2MFILE    Export a PDE from CHEBGUI to a .m file.
%
%   See also: chebgui/export.

% TODO:  Documentation.
% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.


% Extract info from the expInfo struct
deInput = expInfo.deInput;
pdeflag = expInfo.pdeflag;
indVarName = expInfo.indVarName;
sol = expInfo.sol;
sol0 = expInfo.sol0;

% Print description of the problem:
exporter.printDescription(fid, expInfo)


% Print the problem set-up:
exporter.printSetup(fid, expInfo, guifile)

% Print the options set-up:
exporter.printOptions(fid, expInfo)
% 
% % Print the solution step:
% exporter.printSolver(fid, expInfo)
% 



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

% Print the post-solution process:
exporter.printPostSolver(fid, expInfo)

end