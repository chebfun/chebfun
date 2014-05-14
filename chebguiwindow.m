function varargout = chebguiwindow(varargin)
% CHEBGUIWINDOW Driver file for Chebfun's CHEBGUI
%  This m-file populates and controls Chebfun's CHEBGUI.
%
% See also chebgui

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @chebguiwindow_OpeningFcn, ...
    'gui_OutputFcn',  @chebguiwindow_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    warnstate = warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid');
    try
        gui_mainfcn(gui_State, varargin{:});
        warning(warnstate);
    catch ME
        warning(warnstate)
        MEID = ME.identifier;
        if ~isempty(strfind(MEID,'Chebgui:'))
            % These are expected GUI errors. We only show the dialog
            errordlg(cleanErrorMsg(ME.message), 'Chebgui error', 'modal');
            uiwait
            resetComponents(varargin{4});
            
            % If in debug mode, we throw the error to the command window as
            % well
            if get(varargin{4}.menu_debug,'UserData')
                rethrow(ME)
            end
        else
            % Show an error dialog, but also throw the error to the command
            % window
            errordlg(cleanErrorMsg(ME.message), 'Chebgui error', 'modal');
            uiwait
            resetComponents(varargin{4});
            rethrow(ME)
        end
    end
end
% End initialization code - DO NOT EDIT

% --- Executes just before chebguiwindow is made visible.
function chebguiwindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see Output(1-x^2)*exp(-30*(x+.5)^2)Fcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chebguiwindow (see VARARGIN)

% Choose default command line output for chebguiwindow
handles.output = hObject;

chebgui.initialiseFigures(handles)

% Variable that determines whether a solution is available
handles.hasSolution = 0;

% Variable which stores the initial guess/condition
handles.init = [];

% Variable which stores imported variables from workspace
handles.importedVar = struct;


% Get the GUI object from the input argument
if ~isempty(varargin)
    handles.guifile = varargin{1};
else
    cgTemp = chebgui('type','bvp');
    handles.guifile = loadexample(cgTemp,-1); % Load a random example
end
% Create a new structure which contains information about the latest
% solution obtained
handles.latest = struct;

% Store the default length of pausing between plots for BVPs and the
% tolerance in the userdata field of relevant menu objects.
set(handles.menu_odeplottingpause,'UserData','0.5');
set(handles.menu_tolerance,'UserData','1e-10');

% Create UserData for the Fix-Y-axis options (so that we can display
% something if it gets called without selecting a demo).
set(handles.menu_pdefixon,'UserData',{''});

% Populate the Demos menu, but only once (i.e. if user calls chebgui again,
% don't reload the examples).
if isempty(get(handles.menu_demos,'UserData'))
    loaddemo_menu(handles.guifile,handles);
    handles.demosLoaded = 1;
end
% Load the input fields
loadfields(handles.guifile,handles);

% Make sure the GUI starts in the correct mode
switchmode(handles.guifile,handles,handles.guifile.type);

% Get the system font size and store in handles
s = char(com.mathworks.services.FontPrefs.getCodeFont);
if s(end-2) == '='
    fs = round(3/4*str2num(s(end-1)));
else
    fs = round(3/4*str2num(s(end-2:end-1)));
end
set(handles.tempedit,'FontSize',fs);

% set(handles.check_uselatest,'String',{'Use latest';'solution'}); 

% Set the solve button to green
set(handles.button_solve,'String','Solve');
set(handles.button_solve,'BackgroundColor',[43 129 86]/256);

% Ensure that we have a light-grey color in background
set(handles.mainWindow,'BackgroundColor', 0.9*[1 1 1]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chebguiwindow wait for user response (see UIRESUME)
% uiwait(handles.chebguimainwindow);

% --- Outputs from this function are returned to the command line.
function varargout = chebguiwindow_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if nargout > 0,
    varargout{1} = handles.output;
end
% If nargout == 2, return the fll set of handles
if nargout > 1,
    varargout{2} = handles;
end

% -------------------------------------------------------------------------
% ---------- Callback functions for the objects of the GUI ----------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ------------- Functions which call chebgui methods ----------------------
% -------------------------------------------------------------------------

function button_clear_Callback(hObject, eventdata, handles)
if strcmp(get(handles.button_clear,'String'),'Clear all')
    [newGUI handles] = clearGUI(handles.guifile,handles);
    handles.guifile = newGUI;
    guidata(hObject, handles);
elseif strcmp(get(handles.button_clear,'String'),'Pause')
    set(handles.button_clear,'String','Continue');
    set(handles.button_clear,'BackgroundColor',[43 129 86]/256);
    % Re-enable figure buttons
    set(handles.button_figsol,'Enable','on');
    set(handles.button_fignorm,'Enable','on');
else
    % Disable figure buttons
    set(handles.button_figsol,'Enable','off');
    set(handles.button_fignorm,'Enable','off');
    set(handles.button_clear,'String','Pause');
    set(handles.button_clear,'BackgroundColor',[255 179 0]/256);  
end

function button_solve_Callback(hObject, eventdata, handles)
% uicontrol(handles.panel_input)
% figure(handles.chebguimainwindow)
% set(hObject, 'Enable', 'off');
% drawnow;
% set(hObject, 'Enable', 'on');
handles = solveGUI(handles.guifile,handles);
guidata(hObject, handles);

function input_LBC_Callback(hObject, eventdata, handles)
newString = cellstr(get(hObject,'String'));
newString = removeTabs(newString); % Remove tabs
set(hObject,'String',newString);
handles = callbackBCs(handles.guifile,handles,newString,'lbc');
handles.guifile.LBC = newString;
guidata(hObject, handles);


function input_RBC_Callback(hObject, eventdata, handles)
newString = cellstr(get(hObject,'String'));
newString = removeTabs(newString); % Remove tabs
set(hObject,'String',newString);
handles = callbackBCs(handles.guifile,handles,newString,'rbc');
handles.guifile.RBC = newString;
guidata(hObject, handles);

% -------------------------------------------------------------------------
% --------- Functions which do their work without chebgui methods ---------
% -------------------------------------------------------------------------

function dom_left_Callback(hObject, eventdata, handles)
% Store the contents of input1_editText as a string. if the string
% is not a number then input will be empty
input = str2num(get(hObject,'String'));
% Checks to see if input is not numeric or empty. If so, default left end
% of the domain is taken to be -1.
if isempty(input) || isnan(input)
    warndlg('Left endpoint of domain unrecognized, default value -1 used.')
    set(hObject,'String','-1')
end

set(handles.input_GUESS,'Enable','on');
set(handles.toggle_useLatest,'Value',0);
set(handles.toggle_useLatest,'Enable','off');

handles.guifile.DomLeft = get(hObject,'String');
guidata(hObject, handles);


function dom_right_Callback(hObject, eventdata, handles)
input = str2num(get(hObject,'String'));
% Checks to see if input is not numeric or empty. If so, default right end
% of the domain is taken to be 1.
if isempty(input) || isnan(input)
    warndlg('Right endpoint of domain unrecognized, default value 1 used.')
    set(hObject,'String','1')
end

set(handles.input_GUESS,'Enable','on');
set(handles.toggle_useLatest,'Value',0);
set(handles.toggle_useLatest,'Enable','off');

handles.guifile.DomRight = get(hObject,'String');
guidata(hObject, handles);


function input_domain_Callback(hObject, eventdata, handles)
in = get(hObject,'String');
input = str2num(in);
% Checks to see if input is not numeric or empty. If so, default left end
% of the domain is taken to be -1.
if input(1) >= input(end)
    warndlg('Empty domain. Default value [-1,1] used.')
    in = '[-1,1]';
    set(hObject,'String',in);
elseif isempty(input) || any(isnan(input)) || length(input)<2
    warndlg('Domain unrecognized. Default value [-1,1] used.')
    in = '[-1,1]';
    set(hObject,'String',in);
elseif ~any(strfind(in,'['))
    in = ['[' in ']'];
    set(hObject,'String',in);
end

set(handles.input_GUESS,'Enable','on');
set(handles.toggle_useLatest,'Value',0);
set(handles.toggle_useLatest,'Enable','off');

handles.guifile.domain = in;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function input_domain_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% -------------------------------------------------------------------------
% ------- Functions which do their work in a couple of lines of code ------
% -------------------------------------------------------------------------

function input_GUESS_Callback(hObject, eventdata, handles)
newString = cellstr(get(hObject,'String'));

% Remove tabs
newString = removeTabs(newString);
set(hObject,'String',newString);

handles.guifile.init = newString;
if isempty(newString) || (iscell(newString) && numel(newString)==1 && isempty(newString{1}))
    handles.init = '';
    axes(handles.fig_sol);
    cla(handles.fig_sol,'reset');
    guidata(hObject, handles);
    return
end

loadVariables(handles.importedVar)


guidata(hObject, handles);

xtTemp = chebfun('x',str2num(handles.guifile.domain));
% handles.init
if ~exist('r','var'), r = xtTemp; end
if ~exist('x','var'), x = xtTemp; end
if ~exist('t','var'), t = xtTemp; end
% Do something clever with multilines
str = cellstr(get(hObject,'String'));
init = [];
for k = 1:numel(str)
    strk = str{k};
    equalSigns = find(strk=='=');
    if numel(equalSigns) > 1
        error('Chebgui:InitInput','Too many equals signs in input.');
    elseif numel(equalSigns) == 1
        strk = strk(equalSigns+1:end);
    elseif numel(str) > 1
        error('Chebgui:InitInput',...
            'Error constructing initial guess. Input must include the names of the dependent variables, i.e. be on the form "u = %s",...',strk)
    end
        
    strk = deblank(vectorize(strk));
    try
        if ~isempty(strk)
            init = [init eval(strk)];
        end
    catch ME
        error('Chebgui:InitInput',ME.message)
    end
end
handles.init = init;
axes(handles.fig_sol);
plot(handles.init,'linewidth',2)
if ~isempty(handles.guifile.options.fixYaxisLower)
    ylim([str2num(handles.guifile.options.fixYaxisLower) ...
        str2num(handles.guifile.options.fixYaxisUpper)]);
end
if handles.guifile.options.grid, grid on,  end
guidata(hObject, handles);

function loadVariables(importedVar)
fNames = fieldnames(importedVar);
for i=1:length(fNames), assignin('caller',fNames{i},importedVar.(fNames{i})), end


function input_DE_Callback(hObject, eventdata, handles)
str = cellstr(get(hObject,'String'));

% Remove tabs
str = removeTabs(str);

% Update the DE input and store in guifile
set(handles.input_DE,'String',str);
handles.guifile.DE = str;

% Auto PDE and EIG detection
for k = 1:numel(str)
    strk = str{k};
    if any(strfind(strk,'_'))
        if ~get(handles.button_pde,'value')
            handles = switchmode(handles.guifile,handles,'pde');
        end
        break
    elseif any(strfind(strk,'lam') | strfind(strk,'lambda'))
        if ~get(handles.button_eig,'value')
            handles = switchmode(handles.guifile,handles,'eig');
        end
        break
    end
end
guidata(hObject, handles);

function str = removeTabs(str)
for k = 1:numel(str)
    idx = 1;
    strk = str{k};
    while ~isempty(idx)
        idx = strfind(strk,double(9));
        strk(idx) = [];
    end
    str{k} = strk;
end

%% Keypresses
function input_DE_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab'), 
    if strcmp(eventdata.Modifier,'shift')
        if get(handles.button_pde,'value')
            uicontrol(handles.input_timedomain); 
        else
            uicontrol(handles.input_domain); 
        end
    elseif get(handles.button_pde,'value')
        uicontrol(handles.input_LBC); 
    else
        uicontrol(handles.input_BC); 
    end
end
function input_BC_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab')
    if strcmp(eventdata.Modifier,'shift')
        uicontrol(handles.input_DE); 
    else
        uicontrol(handles.input_GUESS); 
    end
end
function input_LBC_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab')
    if strcmp(eventdata.Modifier,'shift')
        uicontrol(handles.input_DE); 
    else
        uicontrol(handles.input_RBC); 
    end
end
function input_RBC_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab')
    if strcmp(eventdata.Modifier,'shift')
        uicontrol(handles.input_LBC); 
    else
        uicontrol(handles.input_GUESS); 
    end
end
function input_GUESS_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab')
    if strcmp(eventdata.Modifier,'shift')
        uicontrol(handles.input_BC); 
    else
        uicontrol(handles.button_solve);
        set(handles.button_solve,'selected','on');
    end
end
function popupmenu_sigma_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab')
    if strcmp(eventdata.Modifier,'shift')
        uicontrol(handles.edit_eigN); 
    else
        uicontrol(handles.button_solve); 
    end
end
function button_solve_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab')
    if strcmp(eventdata.Modifier,'shift')
        if get(handles.button_eig,'value')
            uicontrol(handles.input_BC); 
        else
            uicontrol(handles.input_GUESS); 
        end
    else
        uicontrol(handles.button_clear); 
    end
elseif strcmp(eventdata.Key,'return')
    button_solve_Callback(hObject, eventdata, handles);
end
function button_clear_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab')
    if strcmp(eventdata.Modifier,'shift')
        uicontrol(handles.button_solve); 
    else
        uicontrol(handles.button_export); 
    end
elseif strcmp(eventdata.Key,'return')
    button_clear_Callback(hObject, eventdata, handles);
end
function button_export_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'tab')
    if strcmp(eventdata.Modifier,'shift')
        uicontrol(handles.button_clear); 
    elseif get(handles.button_exportsoln,'enabled')
        uicontrol(handles.button_exportsoln);
    else 
        uicontrol(handles.input_domain);
    end
elseif strcmp(eventdata.Key,'return')
    button_export_Callback(hObject, eventdata, handles);
end
%%

function input_timedomain_Callback(hObject, eventdata, handles)
str = get(hObject,'String');
if iscell(str), str = str{:}; end
num = str2num(str);

options.WindowStyle = 'modal';

% Indicates we had timedomain with negative spacing (0:-.1:1)
while isempty(num) 
    str = inputdlg('Time domain should be a vector of length >2, with positive spacing, at which times solution is returned','Set time domain',...
        1,{'0:.1:1'},options);
    if isempty(str), str = ''; break, end
    str = str{:};
    num = str2num(str);
    
end

while (~isempty(str) && numel(num) < 3)
    h = (num(2)-num(1))/20;
    def = sprintf('%s:%s:%s',num2str(num(1),'%0.0f'),num2str(h,'%0.2g'),num2str(num(2),'%0.0f'));
    str = inputdlg('Time domain should be a vector of length >2 at which times solution is returned','Set time domain',...
        1,{def},options);
    if isempty(str), str = ''; break, end
    str = str{:};
    num = str2num(str);
end
set(handles.input_timedomain,'String',str);   
handles.guifile.timedomain = str;
guidata(hObject, handles);

% -------------------------------------------------------------------------
% -------------------- Unsorted functions  --------------------------------
% -------------------------------------------------------------------------

function button_fignorm_Callback(hObject, eventdata, handles)

% Check the type of the problem
if get(handles.button_ode,'Value');
    latestNorms = handles.latest.norms;
    figure;
    semilogy(latestNorms,'-*','Linewidth',2),title('Norm of updates'), xlabel('Number of iteration')
    if length(latestNorms) > 1
        XTickVec = 1:max(floor(length(latestNorms)/5),1):length(latestNorms);
        set(gca,'XTick', XTickVec), xlim([1 length(latestNorms)]), grid on
    else % Don't display fractions on iteration plots
        set(gca,'XTick', 1)
    end
elseif get(handles.button_pde,'Value');
    u = handles.latest.solution;
    tt = handles.latest.solutionT;
    varnames = handles.varnames;
    
    xLab = handles.indVarName{1};
    tLab = handles.indVarName{2};
    
    if ~iscell(u)
        figure
        waterfall(u,tt,'simple','linewidth',2)
        xlabel(xLab), ylabel(tLab); zlabel(varnames{1});
    else
        figure
        for k = 1:numel(u)
            subplot(1,numel(u),k);
            waterfall(u{k},tt,'simple','linewidth',2)
            xlabel(xLab), ylabel(tLab), zlabel(varnames{k})
            title(varnames{k})
        end
        
        tmp = u{k}(:,1);
        u1 = tmp.vals(1);
        tmp = get(tmp,'vals');
        x1 = tmp(1);
        
        figure
        cols = get(0,'DefaultAxesColorOrder');
        for k = 1:numel(u)
            plot3(x1,tt(1),u1,'linewidth',2,'color',cols(k,:)); hold on
        end
        legend(varnames{:})
        for k = 1:numel(u)
            waterfall(u{k},tt,'simple','linewidth',2,'edgecolor',cols(k,:))
        end
        xlabel(xLab), ylabel(tLab), grid on
    end
else % eigs
    
    figure, h1 = gca;
    
    if strcmp(handles.latest.type,'eig')
        selection = get(handles.iter_list,'Value');
        ploteigenmodes(handles.guifile,handles,selection,[],h1);
    end
    
end

function button_figsol_Callback(hObject, eventdata, handles)
if get(handles.button_ode,'Value');
    latestSolution = handles.latest.solution;
    figure
    plot(latestSolution,'Linewidth',2), title('Solution at end of iteration')
    xlabel(handles.indVarName);
    varnames = handles.varnames;
    if numel(varnames) == 1
        ylabel(varnames)
    end
    % Turn on grid
    if handles.guifile.options.grid, grid on,  end
    if numel(handles.varnames) > 1, legend(handles.varnames), end
elseif get(handles.button_pde,'Value');
    u = handles.latest.solution;
    tt = handles.latest.solutionT;
    
    varnames = handles.varnames;  
    xLab = handles.indVarName{1};
    tLab = handles.indVarName{2};    
    
    titleStr = sprintf('Solution at final time, %s = %f',tLab,tt(end));
    figure
    if ~iscell(u)
        plot(u(:,end),'Linewidth',2)
        xlabel(xLab);
        ylabel(varnames);
        title(titleStr)
    else
        v = chebfun;
        for k = 1:numel(u)
            uk = u{k};
            v(:,k) = uk(:,end);
        end
        plot(v,'Linewidth',2);
        xlabel(xLab);
        legend(varnames);
        title(titleStr)
    end
    % Turn on grid
    if handles.guifile.options.grid, grid on,  end
    
    % Turn on fixed y-limits
    if ~isempty(handles.guifile.options.fixYaxisLower)
        ylim([str2num(handles.guifile.options.fixYaxisLower) ...
            str2num(handles.guifile.options.fixYaxisUpper)]);
    end
else
    figure, h1 = gca;
    if strcmp(handles.latest.type,'eig')
        selection = get(handles.iter_list,'Value');
        ploteigenmodes(handles.guifile,handles,selection,h1,[]);
    end
end


function toggle_useLatest_Callback(hObject, eventdata, handles)
newVal = get(hObject,'Value');

if newVal % User wants to use latest solution
    set(handles.input_GUESS,'String','Using latest solution');
else
    set(handles.input_GUESS,'String','');
    set(handles.input_GUESS,'Enable','On');
    handles.guifile.init = '';
end
guidata(hObject, handles);


% --- Executes when chebguimainwindow is resized.
function chebguimainwindow_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to chebguimainwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function button_ode_Callback(hObject, eventdata, handles)
handles = switchmode(handles.guifile,handles,'bvp');
guidata(hObject, handles);

% --- Executes on button press in button_pde.
function button_pde_Callback(hObject, eventdata, handles)
handles = switchmode(handles.guifile,handles,'pde');
guidata(hObject, handles);

% --- Executes on button press in button_pde.
function button_eig_Callback(hObject, eventdata, handles)
handles = switchmode(handles.guifile,handles,'eig');
guidata(hObject, handles);



% --- Executes on selection change in iter_list.
function iter_list_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns iter_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from iter_list

% Selecting from the list only does something when we are in e-value mode
% Display corresponding e-funcs when clicking on e-value.

% set(handles.fig_sol,'XLimMode','Manual');

if strcmp(handles.latest.type,'eig')
    selection = get(handles.iter_list,'Value');
    ploteigenmodes(handles.guifile,handles,selection);
end

% xlim(handles.fig_norm,xlim_norm); ylim(handles.fig_norm,ylim_norm);
% -------------------------------------------------------------------------
% ---------------------- Other subfunctions -------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% ----------------- All CreateFcn are stored here -------------------------
% -------------------------------------------------------------------------

function dom_left_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dom_right_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function input_DE_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function input_RBC_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function input_GUESS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function input_LBC_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function iter_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fig_sol_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate fig_sol


function fig_norm_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate fig_norm


function input_timedomain_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tempedit_CreateFcn(hObject, eventdata, handles)

function ylim1_CreateFcn(hObject, eventdata, handles)
function ylim1_Callback(hObject, eventdata, handles)
function ylim2_CreateFcn(hObject, eventdata, handles)
function ylim2_Callback(hObject, eventdata, handles)

function button_solve_ButtonDownFcn(hObject, eventdata, handles)
% -------------------------------------------------------------------------
% ----------------- Right-clicking functions ------------------------------
% -------------------------------------------------------------------------

function input_DE_ButtonDownFcn(hObject, eventdata, handles)
chebguiedit('chebguiwindow', handles.chebguimainwindow,'input_DE');
input_DE_Callback(hObject, eventdata, handles); % Go through the callback function
function input_LBC_ButtonDownFcn(hObject, eventdata, handles)
chebguiedit('chebguiwindow', handles.chebguimainwindow,'input_LBC');
input_LBC_Callback(hObject, eventdata, handles); % Go through the callback function
function input_RBC_ButtonDownFcn(hObject, eventdata, handles)
chebguiedit('chebguiwindow', handles.chebguimainwindow,'input_RBC');
input_RBC_Callback(hObject, eventdata, handles); % Go through the callback function
function input_GUESS_ButtonDownFcn(hObject, eventdata, handles)
chebguiedit('chebguiwindow', handles.chebguimainwindow,'input_GUESS');
input_GUESS_Callback(hObject, eventdata, handles); % Go through the callback function

function editfontsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
% ----------------- Callbacks for menu items  ----------------------------
% -------------------------------------------------------------------------

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_opengui_Callback(hObject, eventdata, handles)
[filename pathname filterindex] = uigetfile('*.guifile','Pick a file');
if ~filterindex, return, end
cgTemp = chebgui(fullfile(pathname,filename));
loadfields(cgTemp,handles);
handles.guifile = cgTemp;
if ~isempty(cgTemp.type)
    handles = switchmode(cgTemp,handles,cgTemp.type);
end    
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_savegui_Callback(hObject, eventdata, handles)
[filename pathname filterindex] = uiputfile('*.guifile','Pick a file');
if ~filterindex, return, end

% name = input('What would you like to name this GUI file? ');
name = '';

if get(handles.button_ode,'value')
    demotype = 'bvp';
elseif get(handles.button_pde,'value')
    demotype = 'pde';
else
    demotype = 'eig';
end
dom = get(handles.input_domain,'string');
t = get(handles.input_timedomain,'string');
DE = get(handles.input_DE,'string');
BC = get(handles.input_BC,'string');
LBC = get(handles.input_LBC,'string');
RBC = get(handles.input_RBC,'string');
init = get(handles.input_GUESS,'string');
tol = '';
damping = '';
plotting = '';

DE = mycell2str(DE);
BC = mycell2str(BC);
LBC = mycell2str(LBC);
RBC = mycell2str(RBC);
init = mycell2str(init);

if strcmp(demotype,'pde')
    s = sprintf(['''%s''\n', ...
            'type = ''%s'';\n', ... 
            'domain = ''%s'';\n', ...
            't = ''%s'';\n', ...
            'DE = %s;\n', ...
            'LBC = %s;\n', ...
            'RBC = %s;\n', ...
            'init = %s;\n', ...
            'tol = ''%s'';\n', ...
            'damping = ''%s'';\n', ...
            'plotting = ''%s'';'], ...
    name, demotype, dom, t, DE, LBC, RBC, init, tol, damping, plotting);
else
    s = sprintf(['''%s''\n', ...
            'type = ''%s'';\n', ... 
            'domain = ''%s'';\n', ...
            'DE = %s;\n', ...
            'BC = %s;\n', ...
            'init = %s;\n', ...
            'tol = ''%s'';\n', ...
            'damping = ''%s'';\n', ...
            'plotting = ''%s'';'], ...
    name, demotype, dom, DE, BC, init, tol, damping, plotting);
end
    
fid = fopen(fullfile(pathname,filename),'w+');
fprintf(fid,s);
fclose(fid);

function out = mycell2str(in)
isc = iscell(in);
if ~isc
    in = {in};
    out = '';
else
    out = '{';
end
for k = 1:numel(in)
    if k > 1,  out = [out ' ; ']; end
    out = [out '''' strrep(in{k},'''','''''') ''''];
end
if isc
    out = [out '}'];
end

% --------------------------------------------------------------------
function menu_demos_Callback(hObject, eventdata, handles)


function menu_bvps_Callback(hObject, eventdata, handles)


function menu_ivps_Callback(hObject, eventdata, handles)


function menu_systems_Callback(hObject, eventdata, handles)


function menu_help_Callback(hObject, eventdata, handles)


function menu_openhelp_Callback(hObject, eventdata, handles)
doc('chebgui')


% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_pdesingle_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_pdesystems_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_export_Callback(hObject, eventdata, handles)
% Offer more options if a solution exists.
if handles.hasSolution
    set(handles.menu_exportmatfile,'Enable','on')
    set(handles.menu_exportworkspace,'Enable','on')
else
    set(handles.menu_exportmatfile,'Enable','off')
    set(handles.menu_exportworkspace,'Enable','off')
end

% --------------------------------------------------------------------
function menu_exportmfile_Callback(hObject, eventdata, handles)
export(handles.guifile,handles,'.m')


function menu_exportchebgui_Callback(hObject, eventdata, handles)
export(handles.guifile,handles,'GUI')

function menu_exportworkspace_Callback(hObject, eventdata, handles)
export(handles.guifile,handles,'Workspace')


function menu_exportmatfile_Callback(hObject, eventdata, handles)
export(handles.guifile,handles,'.mat')


function menu_options_Callback(hObject, eventdata, handles)

function menu_close_Callback(hObject, eventdata, handles)
delete(handles.chebguimainwindow)

function mainWindow_DeleteFcn(hObject, eventdata, handles)


function tempedit_Callback(hObject, eventdata, handles)





function menu_odedampednewton_Callback(hObject, eventdata, handles)


function menu_odedampednewtonon_Callback(hObject, eventdata, handles)
handles.guifile.options.damping = '1';
set(handles.menu_odedampednewtonon,'checked','on');
set(handles.menu_odedampednewtonoff,'checked','off');
guidata(hObject, handles);


function menu_odedampednewtonoff_Callback(hObject, eventdata, handles)
handles.guifile.options.damping = '0';
set(handles.menu_odedampednewtonon,'checked','off');
set(handles.menu_odedampednewtonoff,'checked','on');
guidata(hObject, handles);


function menu_odeplotting_Callback(hObject, eventdata, handles)

function menu_pdeplotting_Callback(hObject, eventdata, handles)

function menu_odeplottingon_Callback(hObject, eventdata, handles)
% Obtain length of pause from handles
handles.guifile.options.plotting = get(handles.menu_odeplottingpause,'UserData'); 
set(handles.menu_odeplottingon,'checked','on');
set(handles.menu_odeplottingoff,'checked','off');
guidata(hObject, handles);


function menu_odeplottingoff_Callback(hObject, eventdata, handles)
handles.guifile.options.plotting = 'off'; % Should store value of plotting length
set(handles.menu_odeplottingon,'checked','off');
set(handles.menu_odeplottingoff,'checked','on');
guidata(hObject, handles)


function menu_odeplottingpause_Callback(hObject, eventdata, handles)
options.WindowStyle = 'modal';
valid = 0;
while ~valid
    pauseInput = inputdlg('Length of pause between plots:','Set pause length',...
        1,{get(hObject,'UserData')},options);
    if isempty(pauseInput) % User pressed cancel
        break
    elseif ~isempty(str2num(pauseInput{1})) % Valid input
        valid = 1;
        % Store new value in the UserData of the object
        set(hObject,'UserData',pauseInput{1});
        % Update length of pause in the chebgui object
        handles.guifile.options.plotting = pauseInput{1};
        % Change the marking of the options
        set(handles.menu_odeplottingon,'checked','on');
        set(handles.menu_odeplottingoff,'checked','off');
    else
        f = errordlg('Invalid input.', 'Chebgui error', 'modal');
        uiwait(f); 
    end
end
guidata(hObject, handles)

% --------------------------------------------------------------------
function menu_plotOpt_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_showgrid_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_pdeholdplot_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_pdefix_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_pdeplotfield_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_pdeholdploton_Callback(hObject, eventdata, handles)
handles.guifile.options.pdeholdplot = 1; % Turn holding off
set(handles.menu_pdeholdploton,'checked','on');
set(handles.menu_pdeholdplotoff,'checked','off');
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_pdeholdplotoff_Callback(hObject, eventdata, handles)
handles.guifile.options.pdeholdplot = 0; % Turn holding off
set(handles.menu_pdeholdploton,'checked','off');
set(handles.menu_pdeholdplotoff,'checked','on');
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_pdeplottingon_Callback(hObject, eventdata, handles)
handles.guifile.options.plotting = 'on'; % Obtain length of pause from handles
set(handles.menu_pdeplottingon,'checked','on');
set(handles.menu_pdeplottingoff,'checked','off');
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_pdeplottingoff_Callback(hObject, eventdata, handles)
handles.guifile.options.plotting = 'off'; % Obtain length of pause from handles
set(handles.menu_pdeplottingon,'checked','off');
set(handles.menu_pdeplottingoff,'checked','on');
guidata(hObject, handles);


function menu_tolerance_Callback(hObject, eventdata, handles)
options.WindowStyle = 'modal';
valid = 0;
while ~valid
    tolInput = inputdlg('Tolerance for solution:','Set tolerance',...
        1,{get(hObject,'UserData')},options);
    if isempty(tolInput) % User pressed cancel
        break
    elseif ~isempty(str2num(tolInput{1})) % Valid input
        valid = 1;
        % Store new value in the UserData of the object
        set(hObject,'UserData',tolInput{1});
        % Update the chebgui object
        handles.guifile.tol = tolInput{1};
    else
        f = errordlg('Invalid input.', 'Chebgui error', 'modal');
        uiwait(f); 
    end
end
guidata(hObject, handles)


% --- Executes on button press in check_uselatest.
function check_uselatest_Callback(hObject, eventdata, handles)
% hObject    handle to check_uselatest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_uselatest


% --------------------------------------------------------------------
function menu_showgridon_Callback(hObject, eventdata, handles)
handles.guifile.options.grid = 1; % Turn grid on
set(handles.menu_showgridon,'checked','on');
set(handles.menu_showgridoff,'checked','off');
axes(handles.fig_sol);  grid on
axes(handles.fig_norm); grid on
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_showgridoff_Callback(hObject, eventdata, handles)
handles.guifile.options.grid = 0; % Turn grid off
set(handles.menu_showgridon,'checked','off');
set(handles.menu_showgridoff,'checked','on');
axes(handles.fig_sol); grid off
axes(handles.fig_norm);  grid off
guidata(hObject, handles);


function menu_fixN_Callback(hObject, eventdata, handles)


function menu_annotateplots_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_pdefixon_Callback(hObject, eventdata, handles)
options.WindowStyle = 'modal';
valid = 0;
defaultAnswer = {handles.guifile.options.fixYaxisLower,...
    handles.guifile.options.fixYaxisUpper};
while ~valid
    fixInput = inputdlg({'Lower y-limit:','Upper y-limit:'},'Fix y-axis',...
        1,defaultAnswer,options);
    if isempty(fixInput) % User pressed cancel
        break
    elseif ~isempty(str2num(fixInput{1})) &&  ~isempty(str2num(fixInput{2})) % Valid input
        valid = 1;
        % Store new value in the UserData of the object
        set(hObject,'UserData',fixInput);
        % Update the chebgui object
        handles.guifile.options.fixYaxisLower = fixInput{1};
        handles.guifile.options.fixYaxisUpper = fixInput{2};
    else
        f = errordlg('Invalid input.', 'Chebgui error', 'modal');
        uiwait(f); 
    end
end

% Change checking
set(handles.menu_pdefixon,'checked','on');
set(handles.menu_pdefixoff,'checked','off');
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_pdefixoff_Callback(hObject, eventdata, handles)
% Clear out fix information
handles.guifile.options.fixYaxisLower = '';
handles.guifile.options.fixYaxisUpper = '';

% Change checking
set(handles.menu_pdefixon,'checked','off');
set(handles.menu_pdefixoff,'checked','on');
guidata(hObject, handles);

function sigma_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_eigsscalar_Callback(hObject, eventdata, handles)

% --- Executes on button press in button_export.
function button_export_Callback(hObject, eventdata, handles)
export(handles.guifile,handles,'.m');


% --- Executes on button press in button_realplot.
function button_realplot_Callback(hObject, eventdata, handles)
set(handles.button_realplot,'Value',1)
set(handles.button_imagplot,'Value',0)
set(handles.button_envelope,'Value',0)
selection = get(handles.iter_list,'Value');
ploteigenmodes(handles.guifile,handles,selection)
    
% --- Executes on button press in button_imagplot.
function button_imagplot_Callback(hObject, eventdata, handles)
set(handles.button_realplot,'Value',0)
set(handles.button_imagplot,'Value',1)
set(handles.button_envelope,'Value',0)
selection = get(handles.iter_list,'Value');
ploteigenmodes(handles.guifile,handles,selection)

% --- Executes on button press in button_envelope.
function button_envelope_Callback(hObject, eventdata, handles)
set(handles.button_realplot,'Value',0)
set(handles.button_imagplot,'Value',0)
set(handles.button_envelope,'Value',1)
selection = get(handles.iter_list,'Value');
ploteigenmodes(handles.guifile,handles,selection)


% --------------------------------------------------------------------
function menu_fixNon_Callback(hObject, eventdata, handles)
options.WindowStyle = 'modal';
valid = 0;
defaultAnswer = {handles.guifile.options.fixN};
while ~valid
    fixInput = inputdlg({'Number of gridpoints:'},'Fix space discretisation',...
        1,defaultAnswer,options);
    if isempty(fixInput) % User pressed cancel
        break
    elseif ~mod(str2num(fixInput{1}),1) % Only allow integers 
        valid = 1;
        % Store new value in the UserData of the object
        set(hObject,'UserData',fixInput);
        % Update the chebgui object
        handles.guifile.options.fixN = fixInput{1};
    else
        f = errordlg('Invalid input.', 'Chebgui error', 'modal');
        uiwait(f); 
    end
end

% Change checking
set(handles.menu_fixNon,'checked','on');
set(handles.menu_fixNoff,'checked','off');
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_fixNoff_Callback(hObject, eventdata, handles)
handles.guifile.options.fixN = '';
% Change checking
set(handles.menu_fixNon,'checked','off');
set(handles.menu_fixNoff,'checked','on');
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_sigma_Callback(hObject, eventdata, handles)
% options.WindowStyle = 'modal';
% valid = 0;
% defaultAnswer = {handles.guifile.sigma};
% while ~valid
%     in = inputdlg({'Set sigma, e.g., ''LR'',''SR'',''LM'', ''SM'', etc'},'Set sigma',...
%         1,defaultAnswer,options);
%     if isempty(in) % User pressed cancel
%         break
%     elseif any(strcmpi(in{1},{'LR','SR','LM', 'SM'})) || ~isempty(str2num(in{1}))
%         in = in{1};
%         valid = 1;
%         % Store new value in the UserData of the object
%         set(hObject,'UserData',in);
%         % Update the chebgui object
%         handles.guifile.sigma = in;
%     else
%         f = errordlg('Invalid input.', 'Chebgui error', 'modal');
%         uiwait(f); 
%     end
% end
% guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_eigno_Callback(hObject, eventdata, handles)
% options.WindowStyle = 'modal';
% valid = 0;
% defaultAnswer = {handles.guifile.options.numeigs};
% if isempty(defaultAnswer{:}), defaultAnswer = {'6'}; end % As defined in linop/eigs
% while ~valid
%     in = inputdlg({'Number of eigenvalues to find?'},'Number of eigenvalues?',...
%         1,defaultAnswer,options);
%     if isempty(in) % User pressed cancel
%         break
%     elseif ~mod(str2num(in{1}),1) % Only allow integers
%         in = in{1};
%         valid = 1;
%         % Store new value in the UserData of the object
%         set(hObject,'UserData',in);
%         % Update the chebgui object
%         handles.guifile.options.numeigs = in;
%     else
%         f = errordlg('Invalid input.', 'Chebgui error', 'modal');
%         uiwait(f); 
%     end
% end
% guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_import_Callback(hObject, eventdata, handles)
% Obtain a 'whos' list from the base workspace
variables =  evalin('base','whos');
baseVarNames = {variables.name};
% Create a list dialog
[selection,OK] = listdlg('PromptString','Select variable(s) to import:',...
    'ListString',baseVarNames,'Name','Import variables to GUI',...
    'OKString','Import','ListSize',[160 200],'ListString',baseVarNames);
if ~OK, return, end % User pressed cancel

% Store all the selected variables in the handles
for selCounter = 1:length(selection)
    handles.importedVar.(baseVarNames{selection(selCounter)}) = ...
        evalin('base',baseVarNames{selection(selCounter)});
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_eigssystem_Callback(hObject, eventdata, handles)
% hObject    handle to menu_eigssystem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on chebguimainwindow and none of its controls.
function chebguimainwindow_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to chebguimainwindow (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
edata = eventdata;
if ~isempty(edata.Modifier) && strcmp(edata.Modifier{1},'control')
    switch edata.Key
        case 'return'
            button_solve_Callback(hObject, eventdata, handles);
        case 'e'
            button_export_Callback(hObject, eventdata, handles);
        case {'+','equal'}
            button_fontinc_Callback(hObject, eventdata, handles);
        case {'-','hyphen'}
            button_fontdec_Callback(hObject, eventdata, handles);
        case 'p'
            if any(strcmp(get(handles.button_clear,'String'),{'Pause','Continue'}))
                button_clear_Callback(hObject, eventdata, handles);
            end
    end
end
% PressedKeyNo = double(get(gcbo,'CurrentCharacter'))



function edit_eigN_Callback(hObject, eventdata, handles)
in = get(handles.edit_eigN,'String');
if ~isempty(in) && isempty(str2num(in))
    errordlg('Invalid input. Number of eigenvalues must be an integer.', 'Chebgui error', 'modal');
    set(handles.edit_eigN,'String',handles.guifile.options.numeigs);
else
    handles.guifile.options.numeigs = in;
end
guidata(hObject, handles);

function edit_eigN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit37_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit37_Callback(hObject, eventdata, handles)
in = get(handles.edit_eigN,'String');
if ~isempty(in) && isempty(str2num(in))
    errordlg('Invalid input.', 'Chebgui error', 'modal');
else
    handles.guifile.options.numeigs = in;
end
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_sigma.
function popupmenu_sigma_Callback(hObject, eventdata, handles)
selected = get(hObject,'Value');
switch selected
    case 1
        handles.guifile.sigma = '';
    case 2
        handles.guifile.sigma = 'lm';
    case 3
        handles.guifile.sigma = 'sm';
    case 4
        handles.guifile.sigma = 'lr';
    case 5
        handles.guifile.sigma = 'sr';
    case 6
        handles.guifile.sigma = 'li';
    case 7
        handles.guifile.sigma = 'si';
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_annotateon_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function menu_annotateoff_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function menu_about_Callback(hObject, eventdata, handles)
aboutWindow = dialog('WindowStyle', 'normal', 'Name', 'About chebgui','Position', [500 500 200 200]);
aboutString = sprintf(['Chebgui was developed by Asgeir Birkisson and Nick Hale ',...
    'as an interface to the differential equations in Chebfun.']);
uicontrol(aboutWindow, ... % Text
             'Style','text', ...
             'String', aboutString,...
             'position',[0 -90 200 200])
uicontrol(aboutWindow, ... % Close button
             'Style','pushbutton','String', 'Close','position',[65 10 75 20],'callback',@(a,b,c)delete(aboutWindow))


hPlotAxes = axes(...       % Axes for plotting the selected plot
                 'Parent', aboutWindow, ...
                 'Units', 'normalized', ...
                 'HandleVisibility','callback','position',[0.1 0.6 .8 .25]);

% Plot the logo
f = chebpoly(10);
xx = linspace(-1,.957,1000);
plot(hPlotAxes,xx,f(xx),'linewidth',3)
t = -cos(pi*(2:8)'/10) *0.99;  % cheb extrema (tweaked)
y = 0*t;
h = text(t, y, num2cell(transpose('chebgui')), ...
    'fontsize',16,'hor','cen','vert','mid','parent',hPlotAxes);
flist = listfonts;
k = strmatch('Rockwell',flist);  % 1st choice
k = [k; strmatch('Luxi Serif',flist)];  % 2nd choice
k = [k; strmatch('luxiserif',flist)];  % 2.5th choice
k = [k; strmatch('Times',flist)];  % 3rd choice
if ~isempty(k), set(h,'fontname',flist{k(1)}), end
axis(hPlotAxes,'off')
       
       
function menu_shortcuts_Callback(hObject, eventdata, handles)
% hObject    handle to menu_shortcuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aboutWindow = dialog('WindowStyle', 'normal', 'Name', 'About chebgui','Position', [500 500 200 275]);
aboutString = sprintf(['The following keyboard shortcuts are supported (case insensitive):\n', ...
    'Ctrl + Return = Solve\n', ...
    'Ctrl + P = Pause\n', ...
    'Ctrl + E = Export to .m file\n', ...
    'Ctrl + + = Increase font-size\n', ...
    'Ctrl + - = Decrease font-size',    ]);
uicontrol(aboutWindow, ... % Text
             'Style','text', ...
             'String', aboutString,...
             'position',[0 -40 200 200])
uicontrol(aboutWindow, ... % Close button
             'Style','pushbutton','String', 'Close','position',[65 10 75 20],'callback',@(a,b,c)delete(aboutWindow))


hPlotAxes = axes(...       % Axes for plotting the selected plot
                 'Parent', aboutWindow, ...
                 'Units', 'normalized', ...
                 'HandleVisibility','callback','position',[0.1 0.6 .8 .25]);

% Plot the logo
f = chebpoly(10);
xx = linspace(-1,.957,1000);
plot(hPlotAxes,xx,f(xx),'linewidth',3)
t = -cos(pi*(2:8)'/10) *0.99;  % cheb extrema (tweaked)
y = 0*t;
h = text(t, y, num2cell(transpose('chebgui')), ...
    'fontsize',16,'hor','cen','vert','mid','parent',hPlotAxes);
flist = listfonts;
k = strmatch('Rockwell',flist);  % 1st choice
k = [k; strmatch('Luxi Serif',flist)];  % 2nd choice
k = [k; strmatch('luxiserif',flist)];  % 2.5th choice
k = [k; strmatch('Times',flist)];  % 3rd choice
if ~isempty(k), set(h,'fontname',flist{k(1)}), end
axis(hPlotAxes,'off')


function msg = cleanErrorMsg(msg)
errUsingLoc = strfind(msg,sprintf('Error using'));
if errUsingLoc
    % Trim everything before this
    msg = msg(errUsingLoc:end);
    % Look for first two line breaks (char(10))
    idx = find(msg==char(10),1);
    % Take only what's after the 2nd
    msg = msg(idx+1:end);
end



function input_BC_Callback(hObject, eventdata, handles)
newString = cellstr(get(hObject,'String'));
newString = removeTabs(newString); % Remove tabs
set(hObject,'String',newString);
handles = callbackBCs(handles.guifile,handles,newString,'bc');
handles.guifile.BC = newString;
guidata(hObject, handles);


function input_BC_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_exportsoln.
function button_exportsoln_Callback(hObject, eventdata, handles)
export(handles.guifile,handles,'WorkspaceJustVars')


% --- Executes on button press in togglebutton4.
function togglebutton4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton4


% --- Executes on button press in button_fontinc.
function button_fontinc_Callback(hObject, eventdata, handles)
% hObject    handle to button_fontinc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get a list of all the names in the handles
names = fieldnames(handles);
textLocs = strfind(names,'text_');
inputLocs = strfind(names,'input_');
buttonLocs = strfind(names,'button_');
panelLocs = strfind(names,'panel_');
toggleLocs = strfind(names,'toggle_');

% Deal with variable maximum for input and noninput-type handles
fontsizediff = []; flag = 0;
tmp = strfind(names,'fontsizediff');
if any([tmp{:}])
    fontsizediff = handles.fontsizediff;
end
if isempty(fontsizediff), fontsizediff = 0; end

% Combine all the locations of elements whose font we wish to change
allLocs = strcat(textLocs,inputLocs,buttonLocs,panelLocs,toggleLocs);

for fieldCounter = 1:length(textLocs)
    if ~isempty(allLocs{fieldCounter})
        % Access the field values dynamically using the .( ) call. We
        % access the font size of each element separately as they can have
        % different relative font sizes.
        currFontSize = get(handles.(names{fieldCounter}),'FontSize');
        if currFontSize >= 30
            % Do nothing (global max)
        elseif (isempty(inputLocs{fieldCounter})&&currFontSize >= 16)
            % Size of non input-type handles is further limited
            if ~flag
                % Record this, so a difference doesn't develop
                fontsizediff = min(fontsizediff+1,13); % 13 = 30-16-1
                flag = true; % but only do this once!
            end
        else
            newFontSize = currFontSize + 1;        
            % Update the fontsize, again using the dynamic access
            set(handles.(names{fieldCounter}),'FontSize',newFontSize)     
        end
    end
end
handles.fontsizediff = fontsizediff;  
guidata(hObject, handles);


% --- Executes on button press in button_fontdec.
function button_fontdec_Callback(hObject, eventdata, handles)
% hObject    handle to button_fontdec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get a list of all the names in the handles
names = fieldnames(handles);
textLocs = strfind(names,'text_');
inputLocs = strfind(names,'input_');
buttonLocs = strfind(names,'button_');
panelLocs = strfind(names,'panel_');
toggleLocs = strfind(names,'toggle_');

% Combine all the locations of elements whose font we wish to change
allLocs = strcat(textLocs,inputLocs,buttonLocs,panelLocs,toggleLocs);

% Deal with variable maximum for input and noninput-type handles
fontsizediff = []; flag = 0;
tmp = strfind(names,'fontsizediff');
if any([tmp{:}])
    fontsizediff = handles.fontsizediff;
end
if isempty(fontsizediff), fontsizediff = 0; end

for fieldCounter = 1:length(textLocs)
    if ~isempty(allLocs{fieldCounter})
        % Access the field values dynamically using the .( ) call
        currFontSize = get(handles.(names{fieldCounter}),'FontSize');
        newFontSize = currFontSize;
        if isempty(inputLocs{fieldCounter})
            % Deal with non input-type handles
            if flag
                % We only want update the diff once, so continue 
                continue
            elseif fontsizediff > 0
                % Update the diff and set flag
                fontsizediff = fontsizediff - 1;
                flag = 1;
                continue
            elseif currFontSize > 5
                newFontSize = currFontSize - 1; 
            end
        elseif currFontSize > 5
            % Input-type handles are easier
            newFontSize = currFontSize - 1;   
        end
        % Update the fontsize, again using the dynamic access
        set(handles.(names{fieldCounter}),'FontSize',newFontSize)  
    end
end
handles.fontsizediff = fontsizediff;
guidata(hObject, handles);

function keypress(hObject, eventdata, handles)
newString = cellstr(get(hObject,'String'));
newString = removeTabs(newString); % Remove tabs
set(hObject,'String',newString);
handles = callbackBCs(handles.guifile,handles,newString,'rbc');
handles.guifile.RBC = newString;
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_test_Callback(hObject, eventdata, handles)
T = 0;
folders = {'bvpdemos','pdedemos','eigdemos'};
for k = 1:numel(folders)
    subdir = fullfile(fileparts(which('chebtest')),'@chebgui','private',folders{k});
    subdirlist = dir(subdir);
    subdirnames = { subdirlist.name };
    fprintf([folders{k}(1:3) , '\n']);
    for j = 1:length(subdirnames)
        if subdirlist(j).isdir, continue; end;
        file = fullfile(subdir,subdirnames{j});
        cgTemp = chebgui(file);
        loadfields(cgTemp,handles);
        handles.guifile = cgTemp;
        if ~isempty(cgTemp.type)
            handles = switchmode(cgTemp,handles,cgTemp.type);
        end    
        handles.hasSolution = 0;
        name = subdirnames{j};
        name = strrep(name,'.guifile','');
        fprintf(['  ' , name]);
        try 
            tic
                handles = solveGUI(handles.guifile,handles);
            t = toc;
            guidata(hObject, handles);
            T = T + t;
            if ~handles.hasSolution
                try close('Chebgui error'), end %#ok<TRYNC>
                error('CHEBGUI:test:NoSol','No solution returned');
            end
            fprintf('  passed in %4.4f seconds.\n',t);
        catch ME %#ok<NASGU>
            try close('Chebgui error'), end %#ok<TRYNC>
            fprintf('  FAILED.\n');
        end
    end
end
fprintf('TOTAL TIME = %4.4f.\n',T);


function menu_debug_Callback(hObject, eventdata, handles)
% hObject    handle to menu_debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.menu_debug,'checked'),'on')
    set(handles.menu_debug,'checked','off');
    set(handles.menu_debug,'UserData',0);
else
    set(handles.menu_debug,'checked','on');
    set(handles.menu_debug,'UserData',1);
end
guidata(hObject, handles);


function resetComponents(handles)
% Enable buttons, figures, etc. Set button to 'solve' again
set(handles.button_solve,'String','Solve');
set(handles.button_solve,'BackgroundColor',[43 129 86]/256);
set(handles.button_clear,'String','Clear all');
set(handles.button_clear,'BackgroundColor',get(handles.button_export,'BackgroundColor'));
set(handles.button_figsol,'Enable','on');
set(handles.button_fignorm,'Enable','on');
set(handles.button_exportsoln,'Enable','off');
set(handles.menu_demos,'Enable','on');

% --------------------------------------------------------------------


% --------------------------------------------------------------------


% --- Executes on selection change in popupmenu_bottomFig.
function popupmenu_bottomFig_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_bottomFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_bottomFig contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_bottomFig
newVal = get(hObject,'Value');

% We should not be able to call this method without solution being availble, but
% let's be on the safe side.
if handles.hasSolution
    axes(handles.fig_norm)            
    switch newVal
        case 1
            axes(handles.fig_norm)
            normDelta = handles.normDelta;
            
            % If normDelta is empty, we actually had a linear problem! So don't
            % do anything.
            if ( isempty(normDelta) )
                warndlg(['Problem was linear. Convergence information for' ...
                    ' Newton iteration is not available.'], 'Linear problem');
                set(handles.popupmenu_bottomFig,'Value',2);
                return
            end
            
            semilogy(normDelta,'-*','Linewidth',2),title('Norm of updates'), xlabel('Iteration number')
            if length(normDelta) > 1
                XTickVec = 1:max(floor(length(normDelta)/5),1):length(normDelta);
                set(gca,'XTick', XTickVec), xlim([1 length(normDelta)]), grid on
            else % Don't display fractions on iteration plots
                set(gca,'XTick', 1)
            end
        case 2
            chebpolyplot(handles.latest.solution,'linewidth',2);
    end
end

% --- Executes during object creation, after setting all properties.
function popupmenu_bottomFig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_bottomFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
