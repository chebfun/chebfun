function options = setupODEoptions(guifile, expInfo)
%SETUPODEOPTIONS   Return a CHEBOPPREF specified by CHEBGUI

% Start by constructing a CHEBOPPREF object
options = cheboppref();

% Default tolerance:
defaultTol = options.bvpTol;
tolInput = guifile.tol;
if ( isempty(tolInput) )
    tolNum = defaultTol;
else
    tolNum = str2double(tolInput);
end

% We need a CHEBFUNPREF as well to ensure the tolerance requested is not
% stricter than current CHEBFUN epsilon
chebfunp = chebfunpref;
if ( tolNum < chebfunp.techPrefs.chebfuneps )
    warndlg('Tolerance specified is less than current chebfun epsilon', ...
        'Warning','modal');
    uiwait(gcf)
end

% Set the tolerance for the solution process
options.bvpTol = tolNum;

% Always display iter. information
options.display = 'iter';

% Obtain information about damping and plotting
dampingOnInput = str2num(guifile.options.damping);
plottingOnInput = str2num(guifile.options.plotting);

if ( dampingOnInput )
    options.damping = 1;
else
    options.damping = 0;
end

if ( isempty(plottingOnInput) ) % If empty, we have either 'off' or 'pause'
    if strcmpi(guifile.options.plotting, 'pause')
        options.plotting = 'pause';
    else
        options.plotting = 'off';
    end
else
    options.plotting = plottingOnInput;
end

% Do we want to show grid?
options.grid = guifile.options.grid;

% What discretization do we want?
options.discretization = expInfo.discretization;

% What IVP solver do we want?
options.ivpSolver = guifile.options.ivpSolver;
end
