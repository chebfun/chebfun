function varargout = surf(f, varargin)
%SURF  Plots three cross sections of a CHEBFUN3.
%   SURF(F) or SURF(F, 'SLIDER') plot three cross sections of a CHEBFUN3 
%   object F. It also allows the user to adjust the cross sections using 
%   sliders.
%
%   Use SLICE instead of this, if F is complex-valued.
% 
% See also CHEBFUN3/PLOT, CHEBFUN3/SLICE, CHEBFUN3/ISOSURFACE and 
% CHEBFUN3/SCAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

runSurf3GUI(f);
end % End of function

function runSurf3GUI(f)

dom = f.domain;
numpts = 51;

% Determine proper bounds for axis of figures:
[xx2,yy2,zz2] = ndgrid(linspace(dom(1), dom(2), numpts), ...
    linspace(dom(3), dom(4), numpts), ...
    linspace(dom(5), dom(6), numpts));
F = feval(f, xx2, yy2, zz2);
if ( ~isreal(F) )
    error('CHEBFUN:CHEBFUN3:surf:complex', ...
        'Complex-valued function: Use slice instead of surf.');
    return
end
a = min(F(:))-eps;
b = max(F(:))+eps;

h = instantiateSurf3GUI();
handles = guihandles(h);


xx = dom(1)*ones(numpts);
[yy,zz] = ndgrid(linspace(dom(3), dom(4), numpts), ...
    linspace(dom(5), dom(6), numpts));
v = feval(f,xx,yy,zz);
if ( ~isreal(v) )
    error('CHEBFUN:CHEBFUN3:surf:complex', ...
        'Use slice for complex-valued functions.');
    return
end

set(handles.xSlider, 'Min', dom(1));
set(handles.xSlider, 'Max', dom(2));
set(handles.xSlider, 'Value', xx(1,1));
maxNumber = 15;
set(handles.xSlider, 'SliderStep', [1/maxNumber , 15/maxNumber]);

handles.f = f;
handles.Bnd = [a, b];
handles.yy1 = yy;
handles.zz = zz;

xlim(handles.axes1, [dom(3), dom(4)])
ylim(handles.axes1, [dom(5), dom(6)])
zlim(handles.axes1, [a, b])
axis(handles.axes1,'manual')
surf(handles.axes1, yy, zz, v)
if ( ~isreal(v) )
    error('CHEBFUN:CHEBFUN3:surf:complex', ...
        'Use slice for complex-valued functions.');
    return
end
yy = dom(3)*ones(numpts);
[xx,zz] = ndgrid(linspace(dom(1), dom(2), numpts), ...
    linspace(dom(5), dom(6), numpts));
v = feval(f,xx,yy,zz);
if ( ~isreal(v) )
    error('CHEBFUN:CHEBFUN3:surf:complex', ...
        'Use slice for complex-valued functions.');
    return
end

set(handles.ySlider, 'SliderStep', [1/maxNumber , 15/maxNumber]);
set(handles.ySlider, 'Min', dom(3));
set(handles.ySlider, 'Max', dom(4));
set(handles.ySlider, 'Value', yy(1,1));
handles.xx = xx;

xlim(handles.axes2,[dom(1), dom(2)])
ylim(handles.axes2,[dom(5), dom(6)])
zlim(handles.axes2,[a, b])
axis(handles.axes2,'manual')
surf(handles.axes2, xx, zz, v)
if ( ~isreal(v) )
    error('CHEBFUN:CHEBFUN3:surf:complex', ...
        'Use slice for complex-valued functions.');
    return
end

zz = dom(5)*ones(numpts);
[xx,yy] = ndgrid(linspace(dom(1), dom(2), numpts), ...
    linspace(dom(3), dom(4), numpts));
v = feval(f,xx,yy,zz);
if ( ~isreal(v) )
    error('CHEBFUN:CHEBFUN3:surf:complex', ...
        'Use slice for complex-valued functions.');
    return
end

set(handles.zSlider, 'SliderStep', [1/maxNumber, 15/maxNumber]);
set(handles.zSlider, 'Min', dom(5));
set(handles.zSlider, 'Max', dom(6));
set(handles.zSlider, 'Value', zz(1,1));
handles.yy2 = yy;

xlim(handles.axes3,[dom(1), dom(2)])
ylim(handles.axes3,[dom(3), dom(4)])
zlim(handles.axes3,[a, b])
axis(handles.axes3,'manual')
surf(handles.axes3, xx, yy, v)
if ( ~isreal(v) )
    error('CHEBFUN:CHEBFUN3:surf:complex', ...
        'Use slice for complex-valued functions.');
    return
end

% Update handles structure
guidata(h, handles); 

% Force the figure to clear when another plot is drawn on it so that GUI
% widgets don't linger.  (NB:  This property needs to be reset to 'add' every
% time we change the plot using a slider; otherwise, the slider movement will
% itself clear the figure, which is not what we want.)
set(h, 'NextPlot', 'replacechildren');

end

function h = instantiateSurf3GUI()

% Load up the GUI from the *.fig file.
installDir = chebfunroot();
h = openFigInCurrentFigure([installDir '/@chebfun3/surf.fig']);

% Do any required initialization of the handle graphics objects.
panels = get(h, 'Children'); % 3 panels exist in the surf3.fig.
G1 = get(panels(1), 'Children');
G2 = get(panels(2), 'Children');
G3 = get(panels(3), 'Children');

for i = 1:1:length(G1)
    if ( isa(G1(i), 'matlab.ui.control.UIControl') )
        % Adjust the background colors of the sliders.
        if ( strcmp(G1(i).Style, 'slider') )
            if ( isequal(get(G1(i), 'BackgroundColor'), ...
                    get(0, 'defaultUicontrolBackgroundColor')) )
                set(G1(i), 'BackgroundColor', [.9 .9 .9]);
            end
        end
        G1(i).Callback = @(hObj, data) ...
            zSlider_Callback(hObj, data, guidata(hObj));
    end
end

for i = 1:1:length(G2)
    if ( isa(G2(i), 'matlab.ui.control.UIControl') )
        % Adjust the background colors of the sliders.
        if ( strcmp(G2(i).Style, 'slider') )
            if ( isequal(get(G2(i), 'BackgroundColor'), ...
                    get(0, 'defaultUicontrolBackgroundColor')) )
                set(G2(i), 'BackgroundColor', [.9 .9 .9]);
            end
        end
        G2(i).Callback = @(hObj, data) ...
            ySlider_Callback(hObj, data, guidata(hObj));
    end
end

for i = 1:1:length(G3)
    if ( isa(G3(i), 'matlab.ui.control.UIControl') )
        % Adjust the background colors of the sliders.
        if ( strcmp(G3(i).Style, 'slider') )
            if ( isequal(get(G2(i), 'BackgroundColor'), ...
                    get(0, 'defaultUicontrolBackgroundColor')) )
                set(G3(i), 'BackgroundColor', [.9 .9 .9]);
            end
        end
        G3(i).Callback = @(hObj, data) ...
            xSlider_Callback(hObj, data, guidata(hObj));
    end    
end

% Store handles to GUI objects so that the callbacks can access them. 
guidata(h, guihandles(h));

end

function xSlider_Callback(hObject, eventdata, handles)
% --- Executes on slider movement.
% hObject    handle to xSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nextPlot = get(hObject.Parent.Parent, 'NextPlot');
set(hObject.Parent.Parent, 'NextPlot', 'add');

yy = handles.yy1;
zz = handles.zz;
xslice = get(hObject, 'Value'); %returns position of slider
xx = xslice*ones(size(yy));
v = feval(handles.f, xx, yy, zz);
xlim(handles.axes1, [yy(1,1), yy(end,end)])
ylim(handles.axes1, [zz(1,1), zz(end,end)])
zlim(handles.axes1, [handles.Bnd(1), handles.Bnd(2)])
axis(handles.axes1, 'manual')
surf(handles.axes1, yy, zz, v)

handles.output = hObject;

set(hObject.Parent.Parent, 'NextPlot', nextPlot);

end

function ySlider_Callback(hObject, eventdata, handles)
% --- Executes on slider movement.
% hObject    handle to ySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nextPlot = get(hObject.Parent.Parent, 'NextPlot');
set(hObject.Parent.Parent, 'NextPlot', 'add');

xx = handles.xx;
zz = handles.zz;
yslice = get(hObject, 'Value'); %returns position of slider
yy = yslice*ones(size(xx));
v = feval(handles.f, xx, yy, zz);
xlim(handles.axes2, [xx(1,1), xx(end,end)])
ylim(handles.axes2, [zz(1,1), zz(end,end)])
zlim(handles.axes2, [handles.Bnd(1), handles.Bnd(2)])
axis(handles.axes2, 'manual')
surf(handles.axes2, xx, zz, v)

handles.output = hObject;

set(hObject.Parent.Parent, 'NextPlot', nextPlot);

end

function zSlider_Callback(hObject, eventdata, handles)
% --- Executes on slider movement.
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nextPlot = get(hObject.Parent.Parent, 'NextPlot');
set(hObject.Parent.Parent, 'NextPlot', 'add');

xx = handles.xx;
yy = handles.yy2;
zslice = get(hObject, 'Value'); %returns position of slider
zz = zslice*ones(size(xx));
v = feval(handles.f, xx, yy, zz);
xlim(handles.axes3, [xx(1,1), xx(end,end)])
ylim(handles.axes3, [yy(1,1), yy(end,end)])
zlim(handles.axes3, [handles.Bnd(1), handles.Bnd(2)])
axis(handles.axes3, 'manual')
surf(handles.axes3, xx, yy, v)

handles.output = hObject;

set(hObject.Parent.Parent, 'NextPlot', nextPlot);

end
