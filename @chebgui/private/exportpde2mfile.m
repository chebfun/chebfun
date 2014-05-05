function exportpde2mfile(guifile,pathname,filename)

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

fullFileName = [pathname,filename];
fid = fopen(fullFileName,'wt');

if ispc
    userName = getenv('UserName');
else
    userName = getenv('USER');
end

fprintf(fid,'%% %s - an executable M-file file for solving a PDE.\n',filename);
fprintf(fid,'%% Automatically created from chebfun/chebgui by user %s\n',userName);
fprintf(fid,'%% %s, %s.\n\n',datestr(rem(now,1),13),datestr(floor(now)));

% Extract information from the GUI fields
dom = guifile.domain;
domNum = str2num(dom);
a = num2str(domNum(1)); b = num2str(domNum(end));
deInput = guifile.DE;
lbcInput = guifile.LBC;
rbcInput = guifile.RBC;
deRHSInput = 'u_t';
initInput = guifile.init;
tt = guifile.timedomain;

% Wrap all input strings in a cell (if they're not a cell already)
if isa(deInput,'char'), deInput = cellstr(deInput); end
if isa(lbcInput,'char'), lbcInput = cellstr(lbcInput); end
if isa(rbcInput,'char'), rbcInput = cellstr(rbcInput); end
if isa(deRHSInput,'char'), deRHSInput = cellstr(deRHSInput); end
if isa(initInput,'char'), initInput = cellstr(initInput); end

% deRHSInput = cellstr(repmat('0',numel(deInput),1));
lbcRHSInput = cellstr(repmat('0',numel(lbcInput),1));
rbcRHSInput = cellstr(repmat('0',numel(rbcInput),1));
initRHSInput = cellstr(repmat('0',numel(initInput),1));

% [deString indVarName] = setupFields(deInput,deRHSInput,'DE');
[deString allVarString indVarNameDE pdeVarName pdeflag ignored allVarNames] = setupFields(guifile,deInput,'DE');
if ~any(pdeflag)
    error('Chebgui:chebpde:notapde',['Input does not appear to be a PDE, ', ...
        'or at least is not a supported type.']);
end

% Obtain the independent variable name appearing in the initial condition
if ~isempty(initInput{1})
    [initString ignored indVarNameInit] = setupFields(guifile,initInput,'BC',allVarString);
else
    indVarNameInit = {''};
end

% Make sure we don't have a disrepency in indVarNames. Create a new
% variable indVarName which contains both independent variable, the first
% entry corresponds to space, the second to time.
if ~isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1})
    if strcmp(indVarNameDE{1},indVarNameInit{1})
        indVarName{1} = indVarNameDE{1};
    else
        error('Chebgui:SolveGUIpde','Independent variable names do not agree')
    end
elseif ~isempty(indVarNameInit{1}) && isempty(indVarNameDE{1})
    indVarName{1} = indVarNameInit{1};
elseif isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1})
    indVarName{1} = indVarNameDE{1};
else
    indVarName{1} = 'x'; % Default value
end

if ~isempty(indVarNameDE{2})
    indVarName{2} = indVarNameDE{2}; % This should never be empty for a PDE though.
else
    indVarName{2} = 't'; % Default value
end

if strcmp(indVarName{1},indVarName{2})
     error('Chebgui:SolveGUIpde','The same variable appears to be used as space and time variable');
end

% Create a string with the variables used in the problem
variableString = [',',indVarName{2},',',indVarName{1},','];

xName = indVarName{1};
tName = indVarName{2};

idx = strfind(deString, ')');


% Support for sum and cumsum
sops = {''};
if ~isempty(strfind(deString(idx(1):end),'cumsum('));
    sops = {',sum,cumsum'};
elseif ~isempty(strfind(deString(idx(1):end),'sum('));
    sops = {',sum'};
else
    sops = {''};
end

% Support for periodic boundary conditions
if (~isempty(lbcInput{1}) && strcmpi(lbcInput{1},'periodic')) || ...
        (~isempty(rbcInput{1}) && strcmpi(rbcInput{1},'periodic'))
    lbcInput{1} = []; rbcInput{1} = []; periodic = true;
else
    periodic = false;
end

deString = [deString(1:idx(1)-1), variableString,'diff',sops{:},deString(idx(1):end)];

% Print the PDE
fprintf(fid,'%% Solving\n');
for k = 1:numel(deInput)
    fprintf(fid,'%%   %s,\n',deInput{k});
end
tmpt = eval(tt); 
fprintf(fid,'%% for %s in [%s,%s] and %s in [%s,%s]',xName,a,b,tName,num2str(tmpt(1)),num2str(tmpt(end)));
if ~isempty(lbcInput{1}) || ~isempty(rbcInput{1})
    fprintf(fid,', subject to\n%%');
    if ~isempty(lbcInput{1})
        if numel(lbcInput)==1 && ~any(lbcInput{1}=='=') && ~any(strcmpi(lbcInput{1},{'dirichlet','neumann'}))
            % Sort out when just function values are passed as bcs.
            lbcInput{1} = sprintf('%s = %s',allVarString,lbcInput{1});
        end      
        fprintf(fid,'   ');
        for k = 1:numel(lbcInput)
            fprintf(fid,'%s',lbcInput{k});
            if k~=numel(lbcInput) && numel(lbcInput)>1, fprintf(fid,', '); end
        end
        fprintf(fid,' at %s = % s\n',xName,a);
    end
    if  ~isempty(lbcInput{1}) && ~isempty(rbcInput{1})
        fprintf(fid,'%% and\n%%',xName,a);
    end
    if ~isempty(rbcInput{1})
        if numel(rbcInput)==1 && ~any(rbcInput{1}=='=') && ~any(strcmpi(rbcInput{1},{'dirichlet','neumann'}))
            % Sort out when just function values are passed as bcs.
            rbcInput{1} = sprintf('%s = %s',allVarString,rbcInput{1});
        end
        fprintf(fid,'   ');
        for k = 1:numel(rbcInput)
            fprintf(fid,'%s',rbcInput{k});
            if k~=numel(rbcInput) && numel(rbcInput)>1, fprintf(fid,', '); end
        end
        fprintf(fid,' at %s = %s\n',xName,b);
    end
    fprintf(fid,'\n');
elseif periodic
    fprintf(fid,', subject to periodic boundary conditions.\n\n');
else
    fprintf(fid,'.\n');
end

% fprintf(fid, '%% Create a domain and the linear function on it.\n');
% fprintf(fid,'[d,%s] = domain(%s,%s);\n',indVarName,a,b);
% fprintf(fid,['\n%% Construct a discretisation of the time domain to solve on.\n']);
% fprintf(fid,'t = %s;\n',tt);

fprintf(fid, '%% Create an interval of the space domain,\n');
fprintf(fid,'dom = %s;\n',dom);
fprintf(fid,'%% and a discretisation of the time domain.\n');
fprintf(fid,'%s = %s;\n',tName,tt);

fprintf(fid,'\n%% Make the rhs of the PDE.\n');
fprintf(fid,'pdefun = %s;\n',deString);
if ~all(pdeflag)
    fprintf(fid, ['pdeflag = [',num2str(pdeflag),']; %% Zero when a variable is indep of time.\n']);
end


% Make assignments for left and right BCs.
fprintf(fid,'\n%% Assign boundary conditions.\n');
if ~isempty(lbcInput{1})
    lbcString = setupFields(guifile,lbcInput,'BC',allVarString);
    idx = strfind(lbcString, ')');
    if ~isempty(idx)
        % Support for sum and cumsum
        if ~isempty(strfind(lbcString(idx(1):end),'cumsum('));
            sops = {',sum,cumsum'};
        elseif ~isempty(strfind(lbcString(idx(1):end),'sum('));
            sops = {',sum'};
        else
            sops = {''};
        end
        lbcString = [lbcString(1:idx(1)-1), variableString,'diff', sops{:},lbcString(idx(1):end)];
%             lbcString = strrep(lbcString,'diff','D');
    end
    fprintf(fid,'bc.left = %s;\n',lbcString);
end

if ~isempty(rbcInput{1})
    rbcString = setupFields(guifile,rbcInput,'BC',allVarString);
    idx = strfind(rbcString, ')');
    if ~isempty(idx)
        % Support for sum and cumsum
        if ~isempty(strfind(rbcString(idx(1):end),'cumsum('));
            sops = {',sum,cumsum'};
        elseif ~isempty(strfind(rbcString(idx(1):end),'sum('));
            sops = {',sum'};
        else
            sops = {''};
        end
        rbcString = [rbcString(1:idx(1)-1), variableString,'diff',sops{:},rbcString(idx(1):end)];
%             rbcString = strrep(rbcString,'diff','D');
    end
    fprintf(fid,'bc.right = %s;\n',rbcString);
end

if periodic
    fprintf(fid,'bc = ''periodic'';\n');
end

% Set up the initial condition
fprintf(fid,'\n%% Construct a chebfun of the space variable on the domain,\n');
fprintf(fid,'%s = chebfun(@(%s) %s, dom);\n',xName,xName,xName);
if iscell(initInput) && numel(initInput) > 1
    fprintf(fid,'%% and of the initial conditions.\n');
else
    fprintf(fid,'%% and of the initial condition.\n');
end
if numel(deInput)==1 && ~ischar(deInput)
    % Get the strings of the dependant variable. Just use allVarNames.
    s = allVarNames;
    sol = s{1}; sol0 = [sol '0'];
    findx = strfind(initInput{1},xName);
    initInput = vectorize(char(initInput));
    equalSign = find(initInput=='=',1,'last');
    if ~isempty(equalSign)
        initInput = strtrim(initInput(equalSign+1:end)); 
    end
    if isempty(findx)
        fprintf(fid,'%s = chebfun(%s,dom);\n',sol0,initInput);
    else
        fprintf(fid,'%s = %s;\n',sol0,vectorize(initInput));
    end        
else
    % Get the strings of the dependant variables. Just use allVarNames   
    s = allVarNames;
    % To deal with 'u = ...' etc in intial guesses
    order = []; guesses = []; inits = [];
    % Match LHS of = with variables in allVarNa
    for initCounter = 1:length(initInput)
        currStr = initInput{initCounter};
        equalSign = find(currStr=='=',1,'first');
        currVar = strtrim(currStr(1:equalSign-1));
        match = find(ismember(allVarNames, currVar)==1);
        order = [order;match];
        currInit = strtrim(currStr(1:equalSign-1));
        currGuess = vectorize(strtrim(currStr(equalSign+1:end)));
        guesses = [guesses;{currGuess}];
        inits = [inits;{currInit}];
    end
    [ignored order] = sort(order);
    
    % If the initial guesses are all constants, we need to wrap them in a
    % chebfun call.
    for k = 1:numel(initInput)
        findx = strfind(initInput{k},'x');
        if ~isempty(findx), break, end
    end
    if isempty(findx)
        for k = 1:numel(initInput)
            guesses{k} = sprintf('chebfun(%s,dom)',guesses{k});
        end
    end
    
    % These can be changed
    initText = '0'; sol0 = 'sol0'; sol = 'sol';
    
    for k = 1:numel(initInput)
        fprintf(fid,'%s%s = %s;\n',inits{order(k)},initText,guesses{order(k)});
    end
    fprintf(fid,'%s = [%s%s,',sol0,inits{order(1)},initText);
    for k = 2:numel(initInput)-1
        fprintf(fid,' %s%s,',inits{order(k)},initText);
    end
    fprintf(fid,' %s%s];\n',inits{order(end)},initText);

end

% Option for tolerance
opts = [];
tolInput = guifile.tol;
opts = [opts,'''Eps'',',tolInput];

if ~all(pdeflag)
%     opts = [opts,',''PDEflag'',','[',num2str(pdeflag),']'];
    opts = [opts,',''PDEflag'',','pdeflag'];
end

% Options for plotting
doplot = guifile.options.plotting;
if strcmpi(doplot,'off')
    opts = [opts,',''Plot'',','''off'''];
else
    dohold = guifile.options.pdeholdplot;
    if dohold
        opts = [opts,',''HoldPlot'',','''on'''];
    end
    ylim1 = guifile.options.fixYaxisLower;
    ylim2 = guifile.options.fixYaxisUpper;
    if ~isempty(ylim1) && ~isempty(ylim2)
        opts = [opts,',''Ylim'',[',ylim1,',',ylim2,']'];
    end
%     plotstyle = get(handles.input_plotstyle,'String');
%     if ~isempty(plotstyle)
%         opts = [opts,',''PlotStyle'',''',plotstyle,''''];
%     end
end

% Options for fixed N
if ~isempty(guifile.options.fixN)
    N = str2num(guifile.options.fixN);
    opts = [opts,',''N'',',N];
end        

% Set up preferences
fprintf(fid,'\n%% Setup preferences for solving the problem.\n');
fprintf(fid,'opts = pdeset');
if isempty(opts)
    fprintf(fid,';\n',opts);
else
    fprintf(fid,'(%s);\n',opts);
end

fprintf(fid,['\n%% Solve the problem using pde15s.\n']);
fprintf(fid,'[%s %s] = pde15s(pdefun,%s,%s,bc,opts);\n',indVarName{2},sol,indVarName{2},sol0);

% Conver sol to variable names
if numel(deInput) > 1
    fprintf(fid,'\n%% Recover variable names.\n');
    for k = 1:numel(s)
        fprintf(fid,'%s = %s{%d};\n',s{k},sol,k);
    end
end

% plotting
if numel(deInput) == 1
    fprintf(fid,'\n%% Create plot of the solution.\n');
%     fprintf(fid,'surf(%s,t,''facecolor'',''interp'')\n',sol);
    fprintf(fid,'waterfall(%s,%s,''simple'',''linewidth'',2)\n',sol,indVarName{2});
else
    fprintf(fid,'\n%% Create plots of the solutions.\n');
%     fprintf(fid,'for k = 1:numel(%s)\n',sol);
%     fprintf(fid,'   subplot(1,numel(%s),k)\n',sol);
%      fprintf(fid,'   surf(sol{k},t,''facecolor'',''interp'')\n');
%     fprintf(fid,'end\n');
    M = numel(deInput);
    for k = 1:numel(deInput)
        fprintf(fid,'subplot(1,%d,%d)\n',M,k);
        fprintf(fid,'waterfall(%s,%s,''linewidth'',2)\n',s{k},indVarName{2});
        fprintf(fid,'xlabel(''%s''), ylabel(''%s''), title(''%s'')\n',indVarName{1},indVarName{2},s{k});
    end
end

fclose(fid);
end