function disp(A)
%DISP   Display a CHEBOP by converting it to a pretty-printed string.
%
% See also DISPLAY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

loose = strcmp(get(0, 'FormatSpacing'), 'loose');

if ( isempty(A) )
    fprintf('   (empty chebop)\n');
    if ( loose )
        fprintf('\n')
    end
    
else
    % Store the name of the variable(s) to the operator
    args = [];
    
    % Check whether the chebop is linear in order to be able to display
    % linearity information. 
    if ( isempty(A.op) )
        fprintf('   Empty operator');
    elseif ( islinear(A) )
        fprintf('   Linear operator');
    else
        fprintf('   Nonlinear operator');
    end
    if ( ~isempty(A.op) )
        fprintf(':\n')
        [str, args] = formatOperator(A.op);
        if ( ~isempty(A.opShow) )
            str = A.opShow;
        end
        fprintf('      %s\n', str);
        if ( loose )
            fprintf('\n')
        end
        fprintf('  ');
    end
    fprintf(' operating on chebfun objects defined on:\n')
    Adom = A.domain;
    dom = strtrim(sprintf('%g ', Adom));
    % Replace ' ' with commas in the domain for printing. If the numbers are not
    % all integers, we also want to keep whitespace between the numbers.
    if ( all(floor(Adom) == ceil(Adom)) )
        dom = strrep(dom, ' ', ',');
    else
        dom = strrep(dom, ' ', ', ');
    end
    % Print the domain
    fprintf('      [%s]\n', dom);
    if ( loose )
        fprintf('\n')
    end
    
    if ( ~isempty(A.lbc) || ~isempty(A.rbc) || ~isempty(A.bc) )
        fprintf('   with\n');
        if ( loose )
            fprintf('\n')
        end
    end
    
    if ( ~isempty(A.lbc) )
        fprintf('    left boundary condition(s):\n');
        printLRbc(A.lbc, A.lbcShow, args, loose);
    end
    
    if ( ~isempty(A.rbc) )
        fprintf('    right boundary condition(s):\n');
        printLRbc(A.rbc, A.rbcShow, args, loose);
    end
    
    if ( ~isempty(A.bc) )
        if ( strcmpi(A.bc, 'periodic') )
            fprintf('    periodic boundary conditions.\n');
        else
            fprintf('    constraints:\n');
            if ( ischar(A.lbc) )
                fprintf('      %s\n', A.bc);
            else
                str = stripHandle(func2str(A.bc));
                fprintf('      %s = 0\n', str);
            end
        end
        if ( loose )
            fprintf('\n')
        end
    end
    
    if ( ~isempty(A.maxnorm) )
        mn = A.maxnorm;
        mn = mn(:)'; % Ensure row vector for printing
        if ( length(mn) > 1 )
            maxnormString = ['[' num2str(mn) ']'];
        else
            maxnormString = num2str(mn);
        end
        fprintf('   enforcing a maximum norm of solution: %s\n', maxnormString);
    end

end

end  % End of function DISP().

function str = stripHandle(op)
idx = strfind(op, ')');
str = op(idx(1)+1:end);

end

function printLRbc(bc, bcShow, args, loose)
% Print left or right boundary condition(s)
if ( ischar(bcShow) )
    % We had an assignment such as N.lbc = 'neumann';
    fprintf('      %s\n', bcShow);
elseif ( isnumeric(bcShow) )
    % We had an assignment such as N.lbc = [2; 3];
    if ( ~isempty(args) )
        % How many conditions are we dealing with?
        numBC = length(bcShow);
        
        % Check if we're dealing with a scalar case (e.g. u=2, u'=3) or system
        % case (e.g. u=2, v=3):
        if ( isempty(strfind(args, ',')) )  % Scalar case

            % Extra whitespace to make things align
            if ( numBC < 4)
                extraWS = repmat(' ', 1, numBC - 1);
            else
                extraWS = '    ';
            end
            
            % We always need to print the first condition:
            fprintf('      %s%s = %s\n', args, extraWS, num2str(bcShow(1)));
            % If we got more conditions, print them as well:
            if ( numBC >= 2 )
                fprintf('      %s''%s = %s\n', args, ...
                    extraWS(1:end-1), num2str(bcShow(2)));
            end
            if ( numBC >= 3 )
                fprintf('      %s''''%s = %s\n', args, ...
                    extraWS(1:end-2), num2str(bcShow(3)));
            end
            % Print all remaining conditions:
            for bcCounter = 4:numBC
                fprintf('      %s^(%i) = %s\n', args, ...
                    bcCounter-1, num2str(bcShow(bcCounter)));
            end
        else
           % Vector case. Split args back into individual variables before so we
           % can print them one per line
           args = strsplit(args,',');
           
           for bcCounter = 1:numBC
               fprintf('      %s = %s\n', args{bcCounter}, num2str(bcShow(bcCounter)));
           end
        end
    else
        fprintf('      %s\n', num2str(bcShow));
    end
else
    str = stripHandle(func2str(bc));
    fprintf('      %s = 0\n', str);
end
if ( loose )
    fprintf('\n')
end
end

function [str, args] = formatOperator(op)
% What type of a function did we get passed?
opFunc = functions(op);

% Convert the function to a string:
op = func2str(op);

if ( strcmp(opFunc.type, 'simple') )
    % If the function passed was of the simple type, we simply return the
    % results of func2str above
    str = op;
    
else
    % We got passed an anonymous function. Return a string to print on a nice
    % form.
    idx = strfind(op, ')');
    args = op(3:idx(1)-1);
    % If we have any commas in args, the independent variable (e.g. x) must be
    % included. For nicer output, don't print it
    idxCommas = strfind(args, ',');
    if ( ~isempty(idxCommas) )
        args = args(idxCommas(1)+1:end);
    end
    
    % Isolate the actual operator part from the string
    op = op(idx(1)+1:end);
    
    % If we're dealing with a system, throw away the end [ and ]. Also, split up
    % components into lines.
    if ( op(1) == '[' && op(end) == ']' )
        op = op(2:end-1);
        
        % Replace ; with the ASCII character for new line. Add indentation to
        % make it look nice
        wSpace = repmat(' ', 1, length(args) + 12);
        op = strrep(op, ';', [sprintf('\n'), wSpace]);
    end

    % Output string
    str = [args, ' |--> ', op];
end

end

function s = mat2char(V)
% MAT2CHAR  Convert the varmat to pretty-print string (incl. BCs)

print_bc = true;       % Print BCs or not
defreal = 6;           % Size of matrix to display
dom = get(V,'domain');
numints = numel(dom.endsandbreaks);

if numints == 1 && isempty(V,'bcs')
    print_bc = false;
end

defreal = max(2,ceil(defreal/(numints-1)));
s1 = ['   with n = ',int2str(defreal),' realization:'];

% Evaluate the linop
if print_bc && ~isempty(V,'bcs')
    Vmat = feval(V,defreal,'bc',dom);
else
    Vmat = feval(V,defreal,'nobc',dom);
end

% Some pretty scaling for large/small numbers
M = max(max(abs(Vmat)));
if ( (M > 1e2 || M < 1) && M ~= 0 )
    M = 10^round(log10(M)-1);
    Vmat = Vmat/M;
    s1 = [s1 sprintf('\n    %2.1e * ',M)];
end

if isreal(Vmat)
    s2 = num2str(Vmat,'  %8.4f');
    for k = 1:size(s2,1)
        s2(k,:) = strrep(s2(k,:),' 0.0000','      0');
        if numel(s2(k,:)) > 5,
            s2(k,1:6) = strrep(s2(k,1:6),'0.0000','     0');
        end
        s2(k,:) = strrep(s2(k,:),'-0.0000','      0');
    end
else
    s2 = num2str(Vmat,'%8.2f');
    for k = 1:size(s2,1)
        s2(k,:) = strrep(s2(k,:),'0.00 ','   0 ');
    end
end

% Pad with appropriate spaces
s2 = [ repmat(' ',size(s2,1),5) s2 ];
s = char(s1,'',s2);
end
