function varargout = isosurface(f, varargin)
%ISOSURFACE   Plots an isosurface of a CHEBFUN3.
%   ISOSURFACE(F) plots an isosurface of the CHEBFUN3 object F in a GUI so 
%   that the user can use a slider to change the level.
%
%   ISOSURFACE(F, 'NPTS', N) plots an isosurface of F again with sliders
%   but allows the user to adjust the grid size (the number of points
%   in each direction) used to create the plot.  Default: N = 51.
%
%   ISOSURFACE(F, 'NOSLIDER') plots 3 isosurfaces at three
%   automatically-chosen levels. There is no slider.
%
%   ISOSURFACE(F, LEV) plots the isosurface of F at the level specified by 
%   the scalar LEV. 
%
%   ISOSURFACE(F, LEV, 'NPTS', N) plots the isosurface and also lets
%   the user specify the plotting grid size N.
%
%   ISOSURFACE(F,LEV,...) allows plotting of isosurfaces with specified 
%   level, color and style.
%
%   If F is complex-valued, then its absolute value is plotted.
%
%   Example 1: f = cheb.gallery3('runge');
%            isosurface(f, 0.8)
%
%   Example 2: f = chebfun3(@(x,y,z) tanh(x+y-.3) + cos(x.*y.*z)./(4+x-y-z));
%            isosurface(f)
%
% See also CHEBFUN3/PLOT, CHEBFUN3/SLICE, CHEBFUN3/SCAN and CHEBFUN3/SURF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%   Developer Note: This code might run slowly, not because evaluation is
%   especially slow in Chebfun3, but because isosurface plots can be slow to
%   render in core MATLAB. The time spent on the pure Chebfun side of this 
%   code might be 1% of the time spent just in the line that calls the
%   isosurface command of MATLAB. This is an issue related to how graphics
%   are handled in MATLAB with different graphic cards.

numpts = [];
if ( nargin > 1 )
    for k = 1:nargin-1
        if ( strcmp(varargin{k}, 'npts')) 
            numpts = varargin{k+1};
        end
    end
end
if ( isempty(numpts) )
    numpts = 51;
end

% Call newplot() manually to prepare the figure/axes for plotting because
% built-in isosurface() doesn't.
newplot();

if ( nargin == 1 || (nargin == 3 && strcmp(varargin{1}, 'npts')) )
    % User has specified the size of grid to sample for the isosurface plot.
    runIsosurface3GUI(f, numpts);
else
    holdState = ishold;
    dom = f.domain;
    [xx, yy, zz] = meshgrid(linspace(dom(1), dom(2), numpts), ...
        linspace(dom(3), dom(4), numpts), linspace(dom(5), dom(6), numpts));
    v = feval(f, xx, yy, zz);
    if ( ~isreal(v) )
        v = abs(v);
    end
end

if ( nargin == 2 && strcmp(varargin, 'noslider') ) 
    % Levels are not specified. So, choose 3 levels yourself.
    fMin = min(v(:)); 
    fMax = max(v(:)); 
    fMean = (fMin + fMax)/2;
    isoval1 = (fMin + fMean)/2; 
    isoval2 = fMean; 
    isoval3 = (fMean + fMax)/2;
    p = patch(isosurface(xx, yy, zz, v, isoval1));
    p.FaceColor = 'red'; 
    p.EdgeColor = 'none';
    hold on
    p = patch(isosurface(xx, yy, zz, v, isoval2));
    p.FaceColor = 'green'; 
    p.EdgeColor = 'none';
    p = patch(isosurface(xx, yy, zz, v, isoval3));
    p.FaceColor = 'blue'; 
    p.EdgeColor = 'none';
    
    % Make objects transparent.
    alpha(.4) 
    
    hold off;
    camlight('headlight')
    lighting gouraud
    view(3)
    title('Three level surfaces');
    xlim([dom(1) dom(2)])
    ylim([dom(3) dom(4)])
    zlim([dom(5) dom(6)])    
    legend(sprintf('%3.2f', isoval1), sprintf('%3.2f', isoval2), ...
        sprintf('%3.2f', isoval3));
    axis tight
    if ( ~holdState )
        hold off
    end
     if ( nargout > 0 )
         varargout = {p};
     end    
     return
    
elseif ( nargin == 2 || (nargin == 4 && strcmp(varargin{2}, 'npts')) )
    % Isovalues are given, but colors and style are not specified.
    if iscell(varargin(1))
        isovals = cell2mat(varargin(1));
    else
        isovals = varargin(1);
    end
    
    if ( numel(isovals) == 1 )
         p = patch(isosurface(xx, yy, zz, v, isovals));
         p.FaceColor = 'red'; 
         p.EdgeColor = 'none';
         
         camlight('headlight')
         lighting gouraud
         view(3)
         xlim([dom(1) dom(2)])
         ylim([dom(3) dom(4)])
         zlim([dom(5) dom(6)])
         if ( ~holdState )
             hold off
         end
         if ( nargout > 0 )
             varargout = {p};
         end
    end
    
elseif ( nargin==3 && ~strcmp(varargin{1}, 'npts') ) % Levels, colors and/or style are specified.
        cc = regexp( varargin{2}, '[bgrcmykw]', 'match' );       % color        
        if ( isempty(cc) ) 
            cc{1}= 'g';
        end
    
    if ( numel(varargin(1))==1 ) % Levels
        if iscell(varargin(1))
            isoval1 = cell2mat(varargin(1));
        else
            isoval1 = varargin(1);
        end
    end
    p = patch(isosurface(xx, yy, zz, v, isoval1));
    p.FaceColor = cc{1}; 
    p.EdgeColor = 'none';
    camlight('headlight')
    lighting gouraud
    view(3)
    xlim([dom(1) dom(2)])
    ylim([dom(3) dom(4)])
    zlim([dom(5) dom(6)])
    if ( ~holdState )
        hold off
    end
     if ( nargout > 0 )
         varargout = {p};
     end     
     
elseif ( nargin==5 && strcmp(varargin{3}, 'npts') ) % Levels, colors and/or style are specified.
    cc = regexp( varargin{2}, '[bgrcmykw]', 'match' );       % color
        if ( isempty(cc) ) 
            cc{1}= 'g';
        end
    
    if ( numel(varargin(1))==1 ) % Levels
        if iscell(varargin(1))
            isoval1 = cell2mat(varargin(1));
        else
            isoval1 = varargin(1);
        end
    end
    p = patch(isosurface(xx, yy, zz, v, isoval1));
    p.FaceColor = cc{1}; 
    p.EdgeColor = 'none';
    camlight('headlight')
    lighting gouraud
    view(3)
    xlim([dom(1) dom(2)])
    ylim([dom(3) dom(4)])
    zlim([dom(5) dom(6)])
    axis tight
    if ( ~holdState )
        hold off
    end
     if ( nargout > 0 )
         varargout = {p};
     end
     
end

end % End of function


function runIsosurface3GUI(f, numpts)

h = instantiateIsosurface3();
handles = guihandles(h);

dom = f.domain;
[xx, yy, zz] = meshgrid(linspace(dom(1), dom(2), numpts), ...
    linspace(dom(3), dom(4), numpts), linspace(dom(5), dom(6), numpts));
v = feval(f, xx, yy, zz);
if ( ~isreal(v) )
    v = abs(v);
end
fMin = min(v(:)); 
fMax = max(v(:)); 
isoVal = (fMin + fMax)/2;

set(handles.isosurfaceSlider, 'Min', fMin);
set(handles.isosurfaceSlider, 'Max', fMax);
set(handles.isosurfaceSlider, 'Value', isoVal);

nSteps = 15; % number of steps allowed by the slider
set(handles.isosurfaceSlider, 'SliderStep', [1/nSteps, 1 ]);

% Plot the isosurface
p = patch(isosurface(xx, yy, zz, v, isoVal));
p.FaceColor = 'red';
p.EdgeColor = 'none';
camlight 
lighting gouraud

% Put the current value of the slider on the GUI
set(handles.printedIsoVal, 'String', num2str(isoVal));

view(3)
xlim([dom(1) dom(2)])
ylim([dom(3) dom(4)])
zlim([dom(5) dom(6)])

% Choose default command line output for isosurface
handles.xx = xx;
handles.yy = yy;
handles.zz = zz;
handles.v = v;
handles.isosurfaceSlice = isoVal;
handles.dom = dom;

% Update handles structure
guidata(h, handles);
handles.output = handles.isosurfaceSlider;

end

function h = instantiateIsosurface3()

% Load up the GUI from the *.fig file.
installDir = chebfunroot();
h = openFigInCurrentFigure([installDir '/@chebfun3/isosurface.fig']);

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
            case 'isosurfaceSlider'
                G(i).Callback = @(hObj, data) ...
                    isosurfaceSlider_Callback(hObj, data, guidata(hObj));
        end
    end
end

% Store handles to GUI objects so that the callbacks can access them. 
guidata(h, guihandles(h));

% Force the figure to clear when another plot is drawn on it so that GUI
% widgets don't linger.  (NB:  This property needs to be reset to 'add' every
% time we change the plot using a slider; otherwise, the slider movement will
% itself clear the figure, which is not what we want.)
set(h, 'NextPlot', 'replacechildren');

end

function isosurfaceSlider_Callback(hObject, eventdata, handles)
% --- Executes on zSlider movement.
% hObject    handle to isosurfaceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nextPlot = get(hObject.Parent, 'NextPlot');
set(hObject.Parent, 'NextPlot', 'add');

isoVal = get(hObject, 'Value'); %returns position of slider
dom = handles.dom;

% Clear the old plot
cla(handles.axes1);

p = patch(isosurface(handles.xx, handles.yy, handles.zz, handles.v, isoVal));
p.FaceColor = 'red';
p.EdgeColor = 'none';

camlight 
lighting gouraud

% Put the current value of the slider on the GUI
set(handles.printedIsoVal, 'String', num2str(isoVal));
xlim([dom(1) dom(2)])
ylim([dom(3) dom(4)])
zlim([dom(5) dom(6)])
handles.output = hObject;

set(hObject.Parent, 'NextPlot', nextPlot);

end
