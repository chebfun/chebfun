function varargout = chebguiedit(varargin)
% CHEBGUIEDIT MATLAB code for chebguiedit.fig
%      CHEBGUIEDIT, by itself, creates a new CHEBGUIEDIT or raises the existing
%      singleton*.
%
%      H = CHEBGUIEDIT returns the handle to a new CHEBGUIEDIT or the handle to
%      the existing singleton*.
%
%      CHEBGUIEDIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHEBGUIEDIT.M with the given input arguments.
%
%      CHEBGUIEDIT('Property','Value',...) creates a new CHEBGUIEDIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chebguiedit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chebguiedit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Developers note:
%   A chebguiedit figure gets created when a user right-clicks the input fields
%   of the CHEBGUI figure.

% Edit the above text to modify the response to help chebguiedit

% Last Modified by GUIDE v2.5 31-Jan-2011 09:40:10

%  Copyright 2011 by The University of Oxford and The Chebfun Developers. 
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chebguiedit_OpeningFcn, ...
                   'gui_OutputFcn',  @chebguiedit_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before chebguiedit is made visible.
function chebguiedit_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for chebguiedit

handles.output = hObject;
mainGuiInput = find(strcmp(varargin, 'chebguiwindow'));
if (isempty(mainGuiInput)) ...
    || (length(varargin) <= mainGuiInput) ...
    || (~ishandle(varargin{mainGuiInput+1}))
    disp('-----------------------------------------------------');
    disp('Improper input arguments. ') 
    disp('chebguiedit should only be called from chebgui.')
    disp('-----------------------------------------------------');
end

% Remember the handle, and adjust our position
handles.chebguiwindow = varargin{mainGuiInput+1};

% Obtain handles using GUIDATA with the caller's handle 
mainHandles = guidata(handles.chebguiwindow);
% Set the edit text to the String of the main GUI's button
handles.outputTarget = varargin{3};
set(handles.edit1, 'String', ...
    get(mainHandles.(varargin{3}), 'String'));

% Get the default font size.
set(handles.edit1, 'FontSize', get(mainHandles.tempedit,'FontSize'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chebguiedit wait for user response (see UIRESUME)
uiwait(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = chebguiedit_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} =[];

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)

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

% --- Executes on button press in buttonCancel.
function buttonCancel_Callback(hObject, eventdata, handles)
% Store the used fontsize
updateEditFontSize(hObject, eventdata, handles);
% Resume and close
uiresume(handles.figure1);
delete(handles.figure1)

% --- Executes on button press in buttonClear.
function buttonClear_Callback(hObject, eventdata, handles)
set(handles.edit1,'String','');

function figure_CloseRequestFcn(hObject, eventdata, handles)
% Store the used fontsize
updateEditFontSize(hObject, eventdata, handles);
% Resume and close
uiresume(hObject);

function fontplusbutton_Callback(hObject, eventdata, handles)
% Increase font size.
fs = get(handles.edit1,'FontSize') + 1;
set(handles.edit1,'FontSize',fs);

function fontmbutton_Callback(hObject, eventdata, handles)
% Decrease font size.
fs = get(handles.edit1,'FontSize') - 1;
set(handles.edit1,'FontSize',fs);

function edit1_ButtonDownFcn(hObject, eventdata, handles)
% Is the same as buttonCancel
buttonCancel_Callback(hObject, eventdata, handles)

function updateEditFontSize(hObject, eventdata, handles)
mainHandles = guidata(handles.chebguiwindow);
set(mainHandles.tempedit,'FontSize',get(handles.edit1,'FontSize'));
