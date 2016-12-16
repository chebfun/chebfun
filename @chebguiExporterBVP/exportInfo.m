function expInfo = exportInfo(guifile)
%EXPORTINFO     Extract useful info from a CHEBGUI object for exporting.
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

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract information from the GUI fields
dom = guifile.domain;
deInput = guifile.DE;
bcInput = guifile.BC;
initInput = guifile.init;

% Wrap all input strings in a cell (if they're not a cell already)
if ( isa(deInput, 'char') )
     deInput = cellstr(deInput);
end

if ( isa(bcInput, 'char') )
     bcInput = cellstr(bcInput);
end

if ( isa(initInput, 'char') )
     initInput = cellstr(initInput);
end

% Obtain useful strings describing the differential equation part:
[deString, allVarString, indVarNameDE, dummy, dummy, dummy, allVarNames] = ...
    setupFields(guifile, deInput, 'DE');

% Do some error checking before we do further printing. 

% Ensure that no EIG specific variables appear:
eigNames = {'lambda'; 'lam'; 'l'};
eigMatch = zeros(size(allVarNames));

for eigNameCounter = 1:length(eigNames);
    eigMatch = eigMatch + strcmp(allVarNames, eigNames{eigNameCounter});
end

assert(~any(eigMatch), 'CHEBFUN:CHEBGUIEXPORTER:exportInfo:eig', ...
    ['Problem appears to be an eigenvalue problem. \n Please make sure ' ...
    '''l'', ''lam'' or ''lambda'' do not appear in appear in input.']);

%Check that independent variable name match.

% Obtain the independent variable name appearing in the initial condition:
useLatest = strcmpi(initInput{1}, 'Using latest solution');
if ( ~isempty(initInput{1}) && ~useLatest )
    [dummy, dummy, indVarNameInit] = ...
        setupFields(guifile, initInput, 'INIT', allVarString);
else
    indVarNameInit = {''};
end

% Make sure we don't have a discrepency in indVarNames
if ( ~isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    if ( strcmp(indVarNameDE{1}, indVarNameInit{1}) )
        indVarNameSpace = indVarNameDE{1};
    else
        error('CHEBFUN:CHEBGUIEXPORTERBVP:exportInfo:SolveGUIbvp', 'Independent variable names do not agree')
    end
elseif ( ~isempty(indVarNameInit{1}) && isempty(indVarNameDE{1}) )
    indVarNameSpace = indVarNameInit{1};
elseif ( isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    indVarNameSpace = indVarNameDE{1};
else
    indVarNameSpace = 'x'; % Default value
end

% Replace the 'DUMMYSPACE' variable in the DE field
deString = strrep(deString, 'DUMMYSPACE', indVarNameSpace);
deString = chebguiExporter.prettyPrintFevalString(deString, allVarNames);

% Support for periodic boundary conditions
if ( ~isempty(bcInput{1}) && strcmpi(bcInput{1}, 'periodic') )
    bcInput{1} = [];
    periodic = true;
else
    periodic = false;
end

% What discretization option do we want?
discretization = chebguiExporter.discOption(periodic, dom, ...
    guifile.options.discretization);

% Add spaces to DOM and ALLVARSTRING so it looks nices once we export
dom = strrep(dom, ',', ', ');
allVarString = strrep(allVarString, ',', ', ');

%% Fill up the expInfo struct
expInfo.dom = dom;
expInfo.deInput = deInput;
expInfo.bcInput = bcInput;
expInfo.initInput = initInput;
expInfo.deString = deString;
expInfo.allVarString = allVarString;
expInfo.allVarNames = allVarNames;
expInfo.indVarNameSpace = indVarNameSpace;
expInfo.periodic = periodic;
expInfo.useLatest = useLatest;
expInfo.numVars = length(allVarNames);

% Information related to options set-up
expInfo.tol = guifile.tol;
expInfo.dampingOn = guifile.options.damping;
expInfo.discretization = discretization;
expInfo.plotting = guifile.options.plotting;

end
