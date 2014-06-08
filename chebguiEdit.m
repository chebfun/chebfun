function varargout = chebguiEdit(varargin)
%CHEBGUIEDIT   CHEBGUI edittor.
%   A CHEBGUIEDIT figure gets created when a user right-clicks the input fields
%   of the CHEBGUI figure. It is not intended for use in any other context.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Suppress irritating MLINT warnings: 
%#ok<*INUSL,*DEFNU,*INUSD>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chebguiEdit_OpeningFcn, ...
                   'gui_OutputFcn',  @chebguiEdit_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
if ( nargin && ischar(varargin{1}) )
    gui_State.gui_Callback = str2func(varargin{1});
end

if ( nargout )
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% --- Executes just before chebguiEdit is made visible.
function chebguiEdit_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for chebguiEdit

handles.output = hObject;
mainGuiInput = find(strcmp(varargin, 'chebguiWindow'));
if ( isempty(mainGuiInput) ...
     || (length(varargin) <= mainGuiInput) ...
     || ~ishandle(varargin{mainGuiInput+1}) )
    disp('-----------------------------------------------------');
    disp('Improper input arguments. ') 
    disp('chebguiEdit should only be called from chebgui.')
    disp('-----------------------------------------------------');
end
% TODO: Why not throw an error here?

% Remember the handle, and adjust our position
handles.chebguiwindow = varargin{mainGuiInput+1};

% Obtain handles using GUIDATA with the caller's handle 
mainHandles = guidata(handles.chebguiwindow);
% Set the edit text to the String of the main GUI's button
handles.outputTarget = varargin{3};
set(handles.edit1, 'String', ...
    get(mainHandles.(varargin{3}), 'String'));

% Get the default font size.
set(handles.edit1, 'FontSize', get(mainHandles.tempedit, 'FontSize'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chebguiEdit wait for user response (see UIRESUME)
uiwait(hObject);

end

% --- Outputs from this function are returned to the command line.
function varargout = chebguiEdit_OutputFcn(hObject, eventdata, handles)  
% Get default command line output from handles structure
varargout{1} =[];
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ( ispc && isequal(get(hObject, 'BackgroundColor'), ...
        get(0, 'defaultUicontrolBackgroundColor')) )
    set(hObject, 'BackgroundColor', 'white');
end
end

function edit1_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in buttonOK.
function buttonOK_Callback(hObject, eventdata, handles)
% Obtain handles using GUIDATA with the caller's handle 
mainHandles = guidata(handles.chebguiwindow);
% Set the edit text to the String of the main GUI's button
set(mainHandles.(handles.outputTarget), 'String', ...
    get(handles.edit1, 'String'));
% Store the used fontsize
updateEditFontSize(hObject, eventdata, handles)
% Resume and close
uiresume(handles.figure1);
delete(handles.figure1)
end

% --- Executes on button press in buttonCancel.
function buttonCancel_Callback(hObject, eventdata, handles)
% Store the used fontsize
updateEditFontSize(hObject, eventdata, handles);
% Resume and close
uiresume(handles.figure1);
delete(handles.figure1)
end

% --- Executes on button press in buttonClear.
function buttonClear_Callback(hObject, eventdata, handles) 
set(handles.edit1, 'String', '');
end

function figure_CloseRequestFcn(hObject, eventdata, handles)
% Store the used fontsize
updateEditFontSize(hObject, eventdata, handles);
% Resume and close
uiresume(hObject);
end

function fontplusbutton_Callback(hObject, eventdata, handles)
% Increase font size.
fs = get(handles.edit1, 'FontSize') + 1;
set(handles.edit1, 'FontSize', fs);
end

function fontmbutton_Callback(hObject, eventdata, handles)
% Decrease font size.
fs = get(handles.edit1, 'FontSize') - 1;
set(handles.edit1, 'FontSize', fs);
end

function edit1_ButtonDownFcn(hObject, eventdata, handles)
% Is the same as buttonCancel
buttonCancel_Callback(hObject, eventdata, handles)
end

function updateEditFontSize(hObject, eventdata, handles)
mainHandles = guidata(handles.chebguiwindow);
set(mainHandles.tempedit, 'FontSize', get(handles.edit1, 'FontSize'));
end
