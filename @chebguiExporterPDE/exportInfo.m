function expInfo = exportInfo(guifile)
%EXPORTINFO    Extract useful info from a CHEBGUI object for exporting.
%
% Calling sequence
%
%   EXPINFO = EXPORTINFO(GUIFILE)
%
% where
%
%   GUIFILE:    A CHEBGUI object
%   EXPINFO:    A struct, containing fields with information for exporting to an
%               .m-file.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract information from the GUI fields
dom = guifile.domain;
domNum = str2num(dom);
a = num2str(domNum(1));
b = num2str(domNum(end));
deInput = guifile.DE;
lbcInput = guifile.LBC;
rbcInput = guifile.RBC;
deRHSInput = 'u_t';
initInput = guifile.init;
tt = guifile.timedomain;

% Wrap all input strings in a cell (if they're not a cell already)
if ( isa(deInput, 'char') )
    deInput = cellstr(deInput);
end

if ( isa(lbcInput, 'char') )
    lbcInput = cellstr(lbcInput);
end

if ( isa(rbcInput, 'char') )
    rbcInput = cellstr(rbcInput);
end

if ( isa(deRHSInput, 'char') )
    deRHSInput = cellstr(deRHSInput);
end

if ( isa(initInput, 'char') )
    initInput = cellstr(initInput);
end

lbcRHSInput = cellstr(repmat('0', numel(lbcInput), 1));
rbcRHSInput = cellstr(repmat('0', numel(rbcInput), 1));
initRHSInput = cellstr(repmat('0', numel(initInput), 1));

% Obtain strings for setting up the problem:
[deString, allVarString, indVarNameDE, pdeVarName, pdeflag, ignored, allVarNames] ...
    = setupFields(guifile, deInput, 'DE');
if ( ~any(pdeflag) )
    error('CHEBFUN:CHEBGUIEXPORTERPDE:exportInfo:notPDE', ...
        ['Input does not appear to be a PDE, ', ...
         'or at least is not a supported type.']);
end

% Obtain the independent variable name appearing in the initial condition
if ( ~isempty(initInput{1}) )
    [initString, ignored, indVarNameInit] = ...
        setupFields(guifile, initInput, 'BC', allVarString);
else
    indVarNameInit = {''};
end

% Make sure we don't have a disrepency in indVarNames. Create a new
% variable indVarName which contains both independent variable, the first
% entry corresponds to space, the second to time.
if ( ~isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    if ( strcmp(indVarNameDE{1}, indVarNameInit{1}) )
        indVarName{1} = indVarNameDE{1};
    else
        error('CHEBFUN:CHEBGUIEXPORTERPDE:exportInfo:solveGUIpde', ...
            'Independent variable names do not agree')
    end
elseif ( ~isempty(indVarNameInit{1}) && isempty(indVarNameDE{1}) )
    indVarName{1} = indVarNameInit{1};
elseif ( isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    indVarName{1} = indVarNameDE{1};
else
    indVarName{1} = 'x'; % Default value
end

if ( ~isempty(indVarNameDE{2}) )
    indVarName{2} = indVarNameDE{2}; % This should never be empty for a PDE.
else
    indVarName{2} = 't'; % Default value
end

% Check that we're not using the same variable for space and time!
if ( strcmp(indVarName{1}, indVarName{2}) )
     error('CHEBFUN:CHEBGUIEXPORTERPDE:exportInfo:solveGUIpde', ...
        'The same variable appears to be used as space and time variable');
end

% Ensure that we have a PDE subscript:
if ( ~isempty(indVarNameDE{2}) )
    % Find what the name of the time variable specified is:
    indVarName{2} = indVarNameDE{2};
else
    error('CHEBFUN:CHEBGUIEXPORTERPDE:exportInfo:solveGUIpde', ...
        'No time variable specified, please do so via using subscript.');
end

% Create a string with the variables used in the problem
variableString = [indVarName{2}, ',', indVarName{1}, ','];

% The names of the variables used for space and time:
xName = indVarName{1};
tName = indVarName{2};

% Support for periodic boundary conditions
if ( (~isempty(lbcInput{1}) && strcmpi(lbcInput{1}, 'periodic')) ...
     || (~isempty(rbcInput{1}) && strcmpi(rbcInput{1}, 'periodic')) )
    lbcInput{1} = [];
    rbcInput{1} = [];
    periodic = true;
else
    periodic = false;
end

% Do we have a left BC specified?
if ( ~isempty(lbcInput{1}) )
    lbcString = setupFields(guifile, lbcInput, 'BC', allVarString);
    
    idx = strfind(lbcString, ')');
    if ( ~isempty(idx) )
        % Inject the t and x variables into the anonymous function:
        lbcString = [lbcString(1:2), tName, ',', lbcString(3:end)];

    end

else
    lbcString = [];
end

% Do we have a right BC specified?
if ( ~isempty(rbcInput{1}) )
    rbcString = setupFields(guifile, rbcInput, 'BC', allVarString);

    idx = strfind(rbcString, ')');
    if ( ~isempty(idx) )
        
        % Inject the t and x variables into the anonymous function:
        rbcString = [rbcString(1:2), tName,',', rbcString(3:end)];
        
    end

else
    rbcString = [];
end

% Glue the DESTRING together
deString = ['@(' variableString, deString(3:end)];

% Find what the dependent variables are:
s = allVarNames;
if ( (numel(deInput) == 1) && ~ischar(deInput) )
    
    % Get the strings of the dependent variable. Just use allVarNames.
    sol = s{1};
    sol0 = [sol '0'];
    
else
    
    % These can be changed
    sol0 = 'sol0';
    sol = 'sol';
    
end

%% Fill up the expInfo struct
expInfo.dom = dom;
expInfo.deInput = deInput;

% PDE specific information
expInfo.a = a;
expInfo.b = b;
expInfo.xName = xName;
expInfo.tName = tName;
expInfo.tt = tt;
expInfo.deRHSInput = deRHSInput;
expInfo.lbcRHSInput = lbcRHSInput;
expInfo.rbcRHSInput = rbcRHSInput;
expInfo.initRHSInput = initRHSInput;
expInfo.pdeflag = pdeflag;
expInfo.initInput = initInput;
expInfo.s = s;
expInfo.sol = sol;
expInfo.sol0 = sol0;

% And then some...
expInfo.deString = deString;
expInfo.pdeVarName = pdeVarName;
expInfo.initString = initString;
expInfo.lbcInput = lbcInput;
expInfo.rbcInput = rbcInput;
expInfo.allVarString = allVarString;
expInfo.allVarNames = allVarNames;
expInfo.indVarName = indVarName;
expInfo.periodic = periodic;
expInfo.lbcString = lbcString;
expInfo.rbcString = rbcString;
% Input related to options
expInfo.tol = guifile.tol;
expInfo.doplot = guifile.options.plotting;
expInfo.dohold = guifile.options.pdeholdplot;
expInfo.ylim1 = guifile.options.fixYaxisLower;
expInfo.ylim2 = guifile.options.fixYaxisUpper;
expInfo.fixN = guifile.options.fixN;
expInfo.pdeSolver = guifile.options.pdeSolver;

end
