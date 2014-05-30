function printDescription(fid, expInfo)
% Extract info from the expInfo struct
deInput = expInfo.deInput;
a = expInfo.a;
b = expInfo.b;
xName = expInfo.xName;
tName = expInfo.tName;
tt = expInfo.tt;
lbcInput = expInfo.lbcInput;
rbcInput = expInfo.rbcInput;
allVarString = expInfo.allVarString;
periodic = expInfo.periodic;


% Print the PDE
fprintf(fid, '%% Solving\n');
for k = 1:numel(deInput)
    fprintf(fid, '%%   %s,\n', deInput{k});
end
tmpt = eval(tt);
fprintf(fid, '%% for %s in [%s,%s] and %s in [%s,%s]', xName, a, b, tName, ...
    num2str(tmpt(1)), num2str(tmpt(end)));
if ( ~isempty(lbcInput{1}) || ~isempty(rbcInput{1}) )
    fprintf(fid, ', subject to\n%%');
    if ( ~isempty(lbcInput{1}) )
        if ( (numel(lbcInput) == 1) && ~any(lbcInput{1} == '=') ...
                && ~any(strcmpi(lbcInput{1}, {'dirichlet', 'neumann'})) )
            % Sort out when just function values are passed as bcs.
            lbcInput{1} = sprintf('%s = %s', allVarString, lbcInput{1});
        end
        fprintf(fid, '   ');
        for k = 1:numel(lbcInput)
            fprintf(fid, '%s', lbcInput{k});
            if ( (k ~= numel(lbcInput)) && (numel(lbcInput) > 1) )
                fprintf(fid, ', ');
            end
        end
        fprintf(fid, ' at %s = % s\n', xName, a);
    end
    if  ( ~isempty(lbcInput{1}) && ~isempty(rbcInput{1}) )
        fprintf(fid, '%% and\n%%', xName, a);
    end
    if ( ~isempty(rbcInput{1}) )
        if ( (numel(rbcInput) == 1) && ~any(rbcInput{1} == '=') ...
                && ~any(strcmpi(rbcInput{1}, {'dirichlet', 'neumann'})) )
            % Sort out when just function values are passed as bcs.
            rbcInput{1} = sprintf('%s = %s', allVarString, rbcInput{1});
        end
        fprintf(fid, '   ');
        for k = 1:numel(rbcInput)
            fprintf(fid, '%s', rbcInput{k});
            if ( (k ~= numel(rbcInput)) && (numel(rbcInput) > 1) )
                fprintf(fid,', ');
            end
        end
        fprintf(fid, ' at %s = %s\n', xName, b);
    end
    fprintf(fid, '\n');
elseif ( periodic )
    fprintf(fid, ', subject to periodic boundary conditions.\n\n');
else
    fprintf(fid, '.\n');
end

end