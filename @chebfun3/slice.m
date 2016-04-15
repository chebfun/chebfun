function varargout = slice(f, varargin)
%SLICE  Plots slices of a CHEBFUN3.
%
%   SLICE(F) creates a slice plot of the CHEBFUN3 object F as a GUI so that
%   the user can move sliders.
%
%   SLICE(F, 'NOSLIDER') plots three contours of a CHEBFUN3 object F at 
%   slices corresponding roughly to the location of global maximum of F.
%
%   SLICE(F, sx, sy, sz) allows to set the location of the three slices of 
%   F to be at x = sx, y = sy and z = sz.
% 
%   If f is complex-valued, then phase portraits are plotted.
%
%   Example:
%   slice(f)
%   slice(f, 'noslider')
%   slice(f, -0.7, 0.1, 0.5)
% 
%   See also CHEBFUN3/PLOT, CHEBFUN3/ISOSURFACE, CHEBFUN3/SURF, and SCAN.

if ( nargin == 1 )
        runSlice3GUI(f);
else
    holdState = ishold;
    dom = f.domain;
    numpts = 51;
    [xx, yy, zz] = meshgrid(linspace(dom(1), dom(2), numpts), ...
        linspace(dom(3), dom(4),numpts), linspace(dom(5), dom(6), numpts));
    v = feval(f, xx, yy, zz);
end

if ( nargin == 2 && strcmp(varargin, 'noslider'))
    % Slices are not specified. Choose three yourself.
    %[col,row,tube] = ind2sub(size(v),find(v == max(v(:))));
    % N.B. We have to use meshgrid here not ndgrid. So, rows are the 1st indices.
    if isreal(v)
        [row, col, tube] = ind2sub(size(v), find(v(:) == max(v(:)), 1, 'last'));
        xslice = xx(row, col, tube); 
        yslice = yy(row, col, tube); 
        zslice = zz(row, col, tube); 
        h = slice(xx, yy, zz, v, xslice, yslice, zslice); 
        shading interp
        colorbar
        axis tight
    else
        %% Plot three phase portraits
        % Find the location of maximum of |f|:
        [row, col, tube] = ind2sub(size(v), find(abs(v(:)) == max(abs(v(:))), 1, 'last'));
        xslice = xx(row, col, tube); 
        yslice = yy(row, col, tube); 
        zslice = zz(row, col, tube); 

        h = slice(xx, yy, zz, angle(-v), xslice, yslice, zslice); 
        set(h, 'EdgeColor','none')
        caxis([-pi pi]),
        colormap('hsv')
        axis('equal') 
    end
    str = sprintf('Three slices: x = %2.1g, y = %2.1g, z = %2.1g', ...
        xslice, yslice, zslice);
    title(str)
    if ( ~holdState )
         hold off
     end
     if ( nargout > 0 )
         varargout = {h};
     end
     
elseif ( nargin == 4 )
    % Slices are specified. Get them:
    xslice = cell2mat(varargin(1));
    yslice = cell2mat(varargin(2)); 
    zslice = cell2mat(varargin(3));
    if isreal(v)
        h = slice(xx, yy, zz, v, xslice, yslice, zslice); 
        shading interp
        colorbar, 
        axis tight
    else
        % Plot three phase portraits:
        h = slice(xx, yy, zz, angle(-v), xslice, yslice, zslice); 
        set(h, 'EdgeColor','none')
        caxis([-pi pi]),
        colormap('hsv')
        axis('equal')         
    end    
    if ( ~holdState )
        hold off
    end
     if ( nargout > 0 )
         varargout = {h};
     end
     
end    

end % End of function

function runSlice3GUI(f)

h = instantiateSlice3GUI();
handles = guihandles(h);
    
dom = f.domain;
numpts = 51;

[xx,yy,zz] = meshgrid(linspace(dom(1), dom(2), numpts), ...
    linspace(dom(3), dom(4), numpts), ...
    linspace(dom(5), dom(6), numpts));
v = feval(f, xx, yy, zz);
if isreal(v)
    [row,col,tube] = ind2sub(size(v), find(v(:) == max(v(:)), 1, 'last'));
else
    [row, col, tube] = ind2sub(size(v), find(abs(v(:)) == max(abs(v(:))), 1, 'last'));    
end
xslice = xx(row,col,tube); 
yslice = yy(row,col,tube); 
zslice = zz(row,col,tube); 

set(handles.xSlider, 'Min', dom(1));
set(handles.xSlider, 'Max', dom(2));
set(handles.xSlider, 'Value', xslice);

set(handles.ySlider, 'Min', dom(3));
set(handles.ySlider, 'Max', dom(4));
set(handles.ySlider, 'Value', yslice);

set(handles.zSlider, 'Min', dom(5));
set(handles.zSlider, 'Max', dom(6));
set(handles.zSlider, 'Value', zslice);

nSteps = 15; % number of slices allowed
set(handles.xSlider, 'SliderStep', [1/nSteps , 1 ]);
set(handles.ySlider, 'SliderStep', [1/nSteps , 1 ]);
set(handles.zSlider, 'SliderStep', [1/nSteps , 1 ]);

% Choose default command line output for slice3
handles.xx = xx;
handles.yy = yy;
handles.zz = zz;
handles.xslice = xslice;
handles.yslice = yslice;
handles.zslice = zslice;
handles.v = v;

if isreal(v)
    slice(xx, yy, zz, v, xslice, yslice, zslice)
    shading interp
    %colormap jet
    colorbar
else
    hh = slice(xx, yy, zz, angle(-v), xslice, yslice, zslice); 
    set(hh, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')     
end
% Update handles structure
guidata(h, handles);
handles.output = handles.xSlider;
end

function h = instantiateSlice3GUI()
% Load up the GUI from the *.fig file.
%h = openfig('/Users/user/Desktop/My work/git/chebfun3/@chebfun3/slice.fig', 'invisible');
%h = openfig('/slice.fig', 'invisible');

installDir = chebfunroot();
h = openfig( [installDir '/@chebfun3/slice.fig'], 'invisible');

% Do any required initialization of the handle graphics objects.
G = get(h, 'Children');
for (i = 1:1:length(G))
    if ( isa(G(i), 'matlab.ui.control.UIControl') )
        % Adjust the background colors of the sliders.
        if ( strcmp(G(i).Style, 'slider') )
            if ( isequal(get(G(i), 'BackgroundColor'), ...
                    get(0, 'defaultUicontrolBackgroundColor')) )
                set(G(i), 'BackgroundColor', [.9 .9 .9]);
            end
        end
        % Register callbacks.
        switch ( G(i).Tag )
            case 'xSlider'
                G(i).Callback = @(hObj, data) ...
                    xSlider_Callback(hObj, data, guidata(hObj));
            case 'ySlider'
                G(i).Callback = @(hObj, data) ...
                    ySlider_Callback(hObj, data, guidata(hObj));
            case 'zSlider'
                G(i).Callback = @(hObj, data) ...
                    zSlider_Callback(hObj, data, guidata(hObj));
        end
    end
end

% Add a toolbar to the GUI.
set(h,'toolbar','figure');

% Store handles to GUI objects so that the callbacks can access them. 
guidata(h, guihandles(h));

% Make the GUI window "visible" to the rest of the handle graphics
% system so that things like gcf(), gca(), etc. work properly.
set(h, 'HandleVisibility', 'on');

% Draw the GUI.
set(h, 'Visible', 'on');

end

% --- Executes on xSlider movement.
function xSlider_Callback(hObject, eventdata, handles)
% hObject    handle to xSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xslice = get(hObject, 'Value');         %returns position of slider
yslice = get(handles.ySlider, 'Value'); %returns position of slider
zslice = get(handles.zSlider, 'Value'); %returns position of slider

if ( isreal(handles.v) )
    handles.slice = slice(handles.xx, handles.yy, handles.zz, handles.v, ...
        xslice, yslice, zslice);
    shading interp
    %colormap jet
    colorbar, 
else
    handles.slice = slice(handles.xx, handles.yy, handles.zz, angle(-handles.v), ...
        xslice, yslice, zslice); 
    set(handles.slice, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')
end
handles.output = hObject;

end

% --- Executes on ySlider movement.
function ySlider_Callback(hObject, eventdata, handles)
% hObject    handle to ySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

yslice = get(hObject, 'Value');         %returns position of slider
xslice = get(handles.xSlider, 'Value'); %returns position of slider
zslice = get(handles.zSlider, 'Value'); %returns position of slider

if ( isreal(handles.v) )
    handles.slice = slice(handles.xx, handles.yy, handles.zz, handles.v, ...
        xslice, yslice, zslice);
    shading interp
    %colormap jet
    colorbar, 
else
    handles.slice = slice(handles.xx, handles.yy, handles.zz, angle(-handles.v), ...
        xslice, yslice, zslice); 
    set(handles.slice, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')    
end

handles.output = hObject;

end

% --- Executes on zSlider movement.
function zSlider_Callback(hObject, eventdata, handles)
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

zslice = get(hObject, 'Value');         %returns position of slider
xslice = get(handles.xSlider, 'Value'); %returns position of slider
yslice = get(handles.ySlider, 'Value'); %returns position of slider

if ( isreal(handles.v) )
    handles.slice = slice(handles.xx, handles.yy, handles.zz, handles.v, ...
        xslice, yslice, zslice);
    shading interp
    colorbar
else
    handles.slice = slice(handles.xx, handles.yy, handles.zz, angle(-handles.v), ...
        xslice, yslice, zslice); 
    set(handles.slice, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')    
end

handles.output = hObject;

end