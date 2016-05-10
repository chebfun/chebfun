function explain(f, varargin)
%EXPLAIN    Play a movie explaining the CHEBFUN2 constructor
%
% EXPLAIN(F), where F is a CHEBFUN2, plays a movie about how the CHEBFUN2
% constructor constructed F.
%
% EXPLAIN(F, MODE) speeds up or slow the movie, depending on the value of MODE.
% The options for the value of MODE are as follows:
% 
%   MODE = S, where S is a numerical value, plays the movie S times faster than
%   the regular speed.
%
%   MODE = 'CONTROL' pauses the movie between its sections, waiting for the
%   user to press a key to continue.
%
% Examples:
%   f = cheb.gallery2('smokering'); explain(f)
%   f = chebfun2(@(x,y) franke(x,y)); explain(f, 'control')
%   f = chebfun2(@(x,y) exp(-(x.^2+y.^2)/2)); explain(f, 2)
%   f = chebfun2(@(x,y) cos(x.*y)); explain(f, .5)

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%% Explain only supports chebfun2 on the unit domain
assert(all(f.domain == [-1 1 -1 1]), 'CHEBFUN:CHEBFUN2:explain:domain', ...
    ['Explain only supports chebfun2 objects that are defined on the ' ...
    'unit interval [-1 1 -1 1].']);

%% Create a new figure to plot on:

% Create a new figure, and add space to the bottom for a text box
fig = figure;
shg
% Want thick lines by default
set(fig, 'DefaultLineLineWidth', 3)
figPos = get(fig, 'position');
figPos(3) = figPos(3) + 400;
figPos(4) = figPos(4) + 100;
set(fig, 'position', figPos)
% Put the plotting axes back to normal size (i.e. don't stretch them)
axPos = get(gca, 'position');
axPos(2) = axPos(2) + .2;
axPos(4) = axPos(4) - .2;
set(gca,'position', axPos);
shg

%% Setup the options for playing the movie, in particular, the speed of it.

% Select the mode to play the chebfun2 explanation movie in.
if ( nargin > 1 )
    if ( isnumeric(varargin{1}) )
        mode = 'speed';
        speed = varargin{1};
    else
        mode = varargin{1};
    end
else
    mode = 'speed';
    speed = 1;
end

% Setup some parameters for the screening of the movie. The parameters govern
% the following attributes:
%   PPAN: The speed at which we pan around the function
%   PMOVE: Controls the speed of tilting and drawing of animated lines
%   PFADE: Pause between steps while fading
%   PBREAK:
%   PSECTIONBREAK: The pause between sections of the movie
if ( strcmpi(mode, 'control') )
    ppan   = 0.1;
    pmove  = 0.05;
    pfade  = 0.05;
    pbreak = 0.2;
    psectionbreak = inf;
elseif ( strcmpi(mode, 'speed') )
    c = 1/speed;
    ppan   = 0.1*c;
    pmove  = 0.05*c;
    pfade  = 0.05*c;
    pbreak = 0.2*c;
    psectionbreak = 2*c;
end

% Colors used for plotting lines on the plot. This is useful for old versions of
% Matlab compatability.
co = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];

% For the movie, we only want to look at 9 and 17 point grids (although the
% current preference of MINSAMPLES might be different), as otherwise, it becomes
% too cluttered.
firstSample = 9;
secondSample = 17;

% Store on what rectangle F lives on
rect = f.domain;
%% Print starting text
str = {'This movie explains how a simplified version of the'; 
    'chebfun2 constructor works on your function.'};
textBox = myTextbox(str);
% Scribble to the plot
txt = scribble('Your chebfun2 movie');
plot(txt)
axis([-1 1 -1 1]), axis off
mypause(psectionbreak)

%% Surface plot of the input CHEBFUN2, followed by panning around the function
textBox = myTextbox('This is your chebfun2', textBox);
sf = plot(f);
axis equal, axis off
% Set the correct size of the current axes
set(gca, 'position', axPos);

% Full pan around function:
sf = fullpan(sf, ppan);
% Tilt to a view from above:
sf = tilt(sf, pmove, 1);
% Slowly go from 100% to 50% opaque:
fade(1, .5, sf, pfade);
% Pause before next section:
mypause(psectionbreak)

%% Start explaining the constructor -- Sample on a coarse grid
str = sprintf(['We first sample on a %u x %u Chebyshev grid obtaining ' ...
    'a matrix A.'], firstSample, firstSample);
textBox = myTextbox(str, textBox);
mypause(psectionbreak)

% How many pivots do we want on the coarse grid?
numCoarsePivots = min(2, length(f));

% Do the first ACA
[P, ff, e, scl, sclinf, pts] = myACA(f, firstSample, numCoarsePivots);

%% Plot the first pivot and row/columns on the coarse grid

% Update the textbox
str = {str; 'We then take entry with the largest absolute value.'};
textBox = myTextbox(str, textBox);
mypause(psectionbreak)

% Prepare the figure for plotting
comet = newplot;
% Plot the first pivot and the row/column on top of the grid
val = feval(f, P(1,1), P(1,2));
% Store the pivot location plots
markers = plot3(P(1,1), P(1,2), val, 'Marker', '.', ...
    'Color', co(1,:), 'MarkerSize', 40);

% Update the text box (with pauses in between)
str = 'That''s the big blue dot and the first pivot.';
textBox = myTextbox(str, textBox);
mypause(psectionbreak)
str = {str; 'We take the entry''s value a, column u and row v...'};
textBox = myTextbox(str, textBox);
mypause(psectionbreak)
str = [str; 'and calculate the residual matrix A = A - uv^T/a.'];
textBox = myTextbox(str, textBox);
mypause(2*psectionbreak)

% Draw the first column/row and store the lines
comets = mycomet(P(1,:), comet, co(1,:), pmove/5);
mypause(psectionbreak)
str = sprintf('The norm of the residual is %1.3e.', e(1));
textBox = myTextbox(str, textBox);

% Make an observation if we got passed a function of numerical rank one
if ( length(f) == 1 )
    str = {str; 'Your function is of numerical rank one.';
        'The first step is complete.'};
    textBox = myTextbox(str, textBox);
    mypause(psectionbreak)
end

% Update the plot to show the first residual, rather than original function
delete(sf);
sf = plot(chebfun2(ff{1}, rect));
alpha(sf, .5);
axis([rect, min(min(ff{1})) - .1, scl + .1]);

%% Plot the second pivot on the coarse grid
if ( numCoarsePivots == 2 )
    str = {'We repeat for one more step on the residual matrix.';
        'As before taking the largest absolute value.'};
    textBox = myTextbox(str, textBox);
    
    % Plot the second pivot, update the textbox and pause
    val = feval(f, P(2,1), P(2,2));
    markers = [markers, plot3(P(2,1), P(2,2), val, 'Marker', '.', ...
        'Color', co(2,:), 'MarkerSize', 40)];
    str = [str; 'That''s the red dot.'];
    textBox = myTextbox(str, textBox);
    mypause(psectionbreak)
    
    % Draw the second column/row
    comets = [comets, mycomet(P(2,:), comet, co(2,:), pmove/5)];
    mypause(psectionbreak)
    
    % Update residual information
    str = [str; ...
        sprintf('After the second step the norm of the residual is %1.3e.', ...
        e(2))];
    textBox = myTextbox(str, textBox);
    
    % Make an observation if the function was of numerical rank 2
    if ( length(f) == 2)
        str = [str; 'The first step is complete.'];
        textBox = myTextbox(str, textBox);
        mypause(psectionbreak),
    end
    
    % Update the plot to show the next residual
    delete(sf);
    sf = plot(chebfun2(ff{2}, rect));
    alpha(sf, .5);
    axis([rect, min(min(ff{2})) - .1, scl + .1]);
    mypause(psectionbreak)
end

%% Tidy up after plotting the sampling on coarse grid
delete(comets{1})
if numCoarsePivots > 1
    delete(comets{2})
end
delete(markers), delete(pts);

%% We did not resolve the function with two pivot points, so take more.
if ( length(f) > 2 )
    % Update the text box (with pauses in between)
    str = 'Your function requires more sampling points to resolve.';
    textBox = myTextbox(str, textBox);
    mypause(psectionbreak)
    str = {str; sprintf('So we next sample it on a %u x %u Chebyshev grid.', ...
        secondSample, secondSample)};
    textBox = myTextbox(str, textBox);
    mypause(psectionbreak)
      
    % Update the plot to show the original input chebfun2
    delete(sf)
    sf=plot(f);
    alpha(sf,.5);
    axis equal, axis off
    mypause(psectionbreak)
    
    str = 'We repeat the same process as before, just on a larger matrix...';
    textBox = myTextbox(str, textBox);
    
    % How many pivots do we want on the finer grid?
    numFinePivots = min(6, length(f));
    
    % Do the second ACA
    [P, ff, e, scl, sclinf, pts] = myACA(f, secondSample, numFinePivots);
    
    markers=[];
    comet = newplot;
    comets = [];
    for jj= 1:min(6, length(f))
        val = feval(f,P(jj,1),P(jj,2)) + eps;
        mark = plot3(P(jj,1), P(jj,2), val, 'Marker', '.', ...
            'Color', co(jj,:), 'MarkerSize', 40);
        hold on
        pause(pbreak)
        myline = mycomet(P(jj,:), comet, co(jj,:), pmove/7);
        comets = [myline comets];
        markers = [mark markers];
        delete(sf);
        rect = f.domain;
        sf = plot(chebfun2(ff{jj}, rect)); alpha(sf,.5);
        axis([rect,min(min(ff{jj})) - .2, scl + .1]);
        mypause(psectionbreak)
    end

    for comCounter=1:length(comets)
        delete(comets{comCounter})
    end

    delete(markers)
    delete(pts)
    delete(sf)  % Remove everything so we can pan around.
    sf = plot(f);
    axis([rect, sclinf - .1, scl + .1]);
    axis equal, axis off
else
    delete(sf)  % Remove everything so we can pan around.
    sf = plot(f);
    axis([rect, -scl - .1, scl + .1]);
    axis equal, axis off
end
mypause(psectionbreak)

% If the function is very bad... need to tell people we take the first six
% slices only...
if ( length(f) > 6 )
    str = {str; sprintf(['This process continues and your function was ' ...
        'eventually sampled on a %u x %u Chebyshev grid.'], ...
        length(f.rows(:,1)), length(f.rows(:,1)))};
elseif ( length(f) > 1 )
    str = 'The first stage is done.';
end
textBox = myTextbox(str, textBox);

%% Waterfall plot is what we have drawn.
water = waterfall(f, '-', 'nslices', min(length(f), 6), 'markersize', 35);

% tilt, pan, tilt
sf = tilt(sf, pmove, -1);
sf = fullpan(sf, ppan); pause(pbreak)
sf = tilt(sf, pmove, 1);
campos([0 0 10])
hold on

str = {'To resolve your function we just need to resolve the';
    'selected columns and rows.';
    'Each column and row is a function of just one variable.'};
textBox = myTextbox(str, textBox);
mypause(psectionbreak)

delete(sf);
mypause(psectionbreak)

%% Skeleton version of ACA.
% Get the pivot locations of the actually resolved CHEBFUN2
P = f.pivotLocations;
markers = [];
trun = min(6, length(f));
P = P(1:trun,:);
[xx, yy] = meshgrid(P(:,1), chebpts(33));
mark1 = plot3(xx, yy, scl + 0*xx, '.k', 'MarkerSize', 20);
[xx, yy] = meshgrid(chebpts(33), P(:,2));
mark2 = plot3(xx, yy, scl + 0*xx, '.k', 'MarkerSize', 20);
mypause(psectionbreak)
delete(water)
clf

%% Do coefficient plots of all the columns and rows
str = 'We use Chebfun to approximate it!';
textBox = myTextbox(str, textBox);
mypause(psectionbreak);
str = {str; ['To check the columns and rows of your chebfun2 are resolved ' ...
    'we can look at their Chebyshev coefficients.']};
textBox = myTextbox(str, textBox);
set(gca,'position', axPos);
% Draw coefficient plot on the slices. Here are the columns.
Ccfs = get(f.cols, 'coeffs');
% Only want the first six columns
if ( length(f) > 6 )
    Ccfs(:, 7:length(f)) = [];
end

% Prepare coefficients for plotting:
Ccfs = log10(abs(flipud(Ccfs)));
Ccfs(Ccfs == -inf) = log10(eps);
PivPos = f.pivotLocations;

% First plot columns:
cols = [];
for jj = 1:min(6, length(f))
    xx = PivPos(jj,1)*ones(length(Ccfs), 1);
    yy = linspace(-1, 1, length(Ccfs));
    col=plot3(xx, yy, Ccfs(:,jj), 'Color', co(jj,:));
    hold on
    campos([0 0 10])
    axis equal, axis off
    cols = [col cols];
end

% Now plot rows:
Rcfs = get(f.rows, 'coeffs');

% Only want the first six rows:
if ( length(f) > 6 )
    Rcfs(:, 7:length(f)) = [];
end

Rcfs = log10(abs(Rcfs));
Rcfs(Rcfs == -inf) = log10(eps);

rows=[];
for jj = 1:min(6, length(f))
    xx = linspace(-1, 1, length(Rcfs));
    yy = PivPos(jj,2)*ones(length(Rcfs), 1);
    row = plot3(xx, yy, Rcfs(:,jj), 'Color', co(jj,:));
    rows = [row rows];
end

%% Tilt and camorbit to show the coefficients decaying
% Fix axis zlimit before panning.
zlim([min(min(min(Ccfs)), min(min(Rcfs))) 0]);
axis square;
shg
sf = tilt(sf, pmove, -1);
pause(pbreak);

k = 36;
for i = 1:k
    camorbit(1, 0, 'data', [0 0 1])
    camorbit(0, 1 ,'data', [0 1 0])
    pause(pmove)
    drawnow
end

str = 'This shows that the rows are resolved.';
textBox = myTextbox(str, textBox);
mypause(psectionbreak)

for i = 1:k
    camorbit(-1, 0, 'data', [0 0 1])
    camorbit(0, -1, 'data', [0 1 0])
    pause(pmove)
    drawnow
end

pause(pbreak)

for i = 1:k
    camorbit(-1.5, 0, 'data', [0 0 1])
    camorbit(0, 1, 'data', [0 1 0])
    pause(pmove)
    drawnow
end

str = 'This shows that the columns are resolved.';
textBox = myTextbox(str, textBox);
mypause(psectionbreak)

for i = 1:k
    camorbit(1.5, 0, 'data', [0 0 1])
    camorbit(0, -1, 'data', [0 1 0])
    pause(pmove)
    drawnow
end

pause(pbreak)
sf = tilt(sf, pmove, 1);


% Put back the pivot locations.
P = PivPos;
U = f.pivotValues;
markers = [];
for jj = 1:min(6, length(f))
    val = feval(f, P(jj,1), P(jj,2)) + eps;
    mark = plot3(P(jj,1), P(jj,2), val, 'Color', co(jj,:), 'Marker', '.', ...
        'MarkerSize', 40);
    markers = [mark markers];
end
campos([0 0 10]);

mypause(psectionbreak)
clf

%% Show how the chebfun2 is stored
txt = scribble('How is it stored?');
plot(txt);
axis([-1 1 -1 1])
axis off
mypause(psectionbreak)
clf

% Make a plot that looks a bit like plot(f, '.-')
crosses = PivPos(1:min(6, length(f)),:);
cols = line([crosses(:,1) crosses(:,1)].', [-1 1]);
hold on
rows = line([-1 1], [crosses(:,2),crosses(:,2)].');
axis(2.3*[-1 1 -1 1]), axis square, axis off

% Again plot pivot positions.
P = PivPos;
U = f.pivotValues;
markers = [];
for jj = 1:min(6, length(f))
    mark = plot(P(jj,1), P(jj,2), 'Marker', '.', 'MarkerSize', 40, ...
        'Color', co(jj,:));
    markers = [mark markers];
end

% Slowly move columns and rows to low rank storage structure.
if length(f) > 1
    cc = linspace(-1.75, -1, min(length(f), 6));
    rr = linspace(.25, 1, min(length(f), 6));
    rr = rr(end:-1:1);
    dd = linspace(-.6, .15, min(length(f), 6));
else
    cc = -1;
    rr = 1;
    dd = -.6;
end

str = 'It''s stored as two chebfun quasi-matrices just like this:';
textBox = myTextbox(str, textBox);
pause(pbreak)

for t = 0:.1:1
    delete(cols)
    delete(rows)
    delete(markers)

    cols = [];
    for jj= 1:min(6, length(f))
        st = t*(cc(jj) - crosses(jj, 1));
        col = line([crosses(jj,1) + st, crosses(jj,1) + st].', [-1 1], ...
            'Color', co(jj,:));
        hold on
        cols = [col, cols];
    end

    rows = [];
    for jj= 1:min(6, length(f))
        st = t*(rr(jj) - crosses(jj,2));
        row = line([-1 + 1.5*t, 1 + 1.5*t], ...
            [crosses(jj,2) + st, crosses(jj, 2) + st].', ...
            'Color', co(jj,:));
        rows = [row, rows];
    end

    markers = [];
    for jj= 1:min(6, length(f))
        stx = t*(dd(jj) - P(jj,1)); sty = t*(rr(jj) - P(jj,2));
        mark = plot(P(jj,1) + stx, P(jj,2) + sty, 'Marker', '.', ...
            'MarkerSize', 40, 'Color', co(jj,:));
        markers = [mark markers];
    end

    pause(pbreak);
end

% Create textbox for equation:
delete(textBox);
annotation(gcf, 'textbox', ...
    [0.243441860465116 0.701923076923077 0.764988372093023 0.0913461538461542], ...
    'String', {'f(x,y) = C(y)         U          R(x)'}, ...
    'FontSize', 20, ...
    'FitBoxToText', 'off', ...
    'LineStyle', 'none', ...
    'LineWidth', 3);

str = 'The end. Thanks for watching!';
pause(pbreak)
textBox = myTextbox(str, textBox);

end


function h = fullpan(h, p)
% Do a full pan

% Set the camera view angle to manual mode, so that the size of the plot
% doesn't keep changing as we rotate
set(gca, 'cameraviewanglemode', 'manual')
axis vis3d
k = 36;
for i = 1:k
    % Spin the camera around the z-axis (the third and fourth arguments to
    % camorbit specify what axis we want to spin around)
    camorbit(10, 0, 'data', [0 0 1])
    pause(p)
    drawnow
end

end

function h = tilt(h, p, whichway)
%TILT
set(gca, 'cameraviewanglemode', 'auto')
k = 36;
s = whichway;
for i = 1:k
    camorbit(0, s*1.665, 'data', [0 1 1])
    camorbit(s*1.04, 0, 'data', [0 0 1])
    pause(p)
    drawnow
end

end

function fade(alpha1, alpha2, h, p)
%FADE
for a = alpha1:-.1:alpha2
    alpha(h, a);
    pause(p);
end

end

function h = mycomet(pivot, h, color, p)
% Draw lines, a row and a column out from a pivot
t = 0:.01:1;
x = zeros(1, length(t));
xeven = pivot(1) + t(1:2:end)*(1 - pivot(1));
xodd = pivot(1) + t(2:2:end)*(-1 - pivot(1));
x(1:2:end) = xeven;
x(2:2:end) = xodd;
y = pivot(2)*ones(length(t), 1);

if ( verLessThan('matlab', '8.4') )
    body = line('parent', h, 'color', color, 'linestyle', '-', ...
        'LineWidth', 3, 'erase', 'none', 'xdata', [], 'ydata', []);

    for jj = 1
        set(body, 'xdata', x(jj), 'ydata', y(jj)),
        drawnow
        hold on
    end

    for jj = 2:length(x)
        j = jj-1:jj;
        set(body, 'xdata', x(j), 'ydata', y(j))
        drawnow
        pause(p)
    end
    
    body2 = line('parent', h, 'color', color, 'linestyle', '-', ...
        'LineWidth', 2, 'erase','none', 'xdata', [], 'ydata', []);
    x = pivot(1)*ones(length(t), 1);
    yeven = pivot(2) + t(1:2:end)*(1 - pivot(2));
    yodd = pivot(2) + t(2:2:end)*(-1 - pivot(2));
    y(1:2:end) = yeven;
    y(2:2:end) = yodd;

    for jj = 1
        set(body2, 'xdata', x(jj), 'ydata', y(jj))
        drawnow
    end

    for jj = 2:length(x)
        j = jj-1:jj;
        set(body2, 'xdata', x(j), 'ydata', y(j))
        drawnow
        pause(p)
    end
else
    body = animatedline('parent', h, 'color', color, 'linestyle', '-', ...
        'LineWidth', 3);
    body2 = animatedline('parent', h, 'color', color, 'linestyle', '-', ...
        'LineWidth', 3);
    for k = 1:length(x);
        addpoints(body, x(k), y(k))
        pause(p)
        drawnow update
    end
    
    x = pivot(1)*ones(length(t), 1);
    yeven = pivot(2) + t(1:2:end)*(1 - pivot(2));
    yodd = pivot(2) + t(2:2:end)*(-1 - pivot(2));
    y(1:2:end) = yeven;
    y(2:2:end) = yodd;

    for k = 1:length(x)
        addpoints(body2, x(k), y(k))
        pause(p)
        drawnow update
    end
end


h = {[body body2]};

end


function mypause(n)
%MYPAUSE    Pause the chebfun2 explanation movie
%
% Input:
%   n: If n < Inf, the movie gets paused for n seconds. If n == Inf, the movie
%   gets paused until the user presses a key.

% Pause depending on the input.
if ( n == inf ) 
    pause
else
    pause(n)
end

end

function textBox = myTextbox(str, oldBox)
%MYTEXTBOX    Draw a textbox to the figure.
%
% Inputs:
%   STR: String to be displayed within the textbox.
%   OLDBOX: A MATLAB graphics object of the currently displayed textbox.
% Output:
%   TEXTBOX: A MATLAB graphics object, the next textbox.

% If OLDBOX got passed in, delete it.
if (nargin > 1 )
    try
        delete(oldBox)
    end
end

% Draw a new textbox.
textBox = annotation('textbox', [0.05, 0.025, 0.9, 0.25], ...
           'String', str, 'fontsize', 18, 'linewidth', 1);

end

function [P, ff, e, scl, sclinf, pts] = myACA(f, nGridPoints, nIter)
% MYACA    A basic ACA algorithm for approximating a 2D function

% Sample the function on a finer Chebyshev grid so that we can compute (and
% plot) the residual after taking one step.
nPlotGrid = 200;
[xx, yy] = meshgrid(chebpts(nPlotGrid));
B = feval(f, xx, yy);

% We first sample on a 9x9 grid.
[xx, yy] = meshgrid(chebpts(nGridPoints));
A = feval(f, xx, yy);
scl = max(max(A));
sclinf = min(min(A));

% Find out what z location we need to plot the grid
zz = (max(max(A)) - 1e-1)*ones(nGridPoints^2, 1);
hold on
pts = plot3(xx(:), yy(:), zz, 'Marker', '.', 'MarkerSize', 20, ...
    'Color', 'k', 'LineStyle', 'none');

% Find maximum on the grid, calculate the residual and then take next maximum.
xpts = chebpts(nGridPoints);

% Store the first one or two pivot locations, residual and error after each step:
P = zeros(nIter, 2);
ff = cell(nIter, 1);
e = zeros(nIter, 1);
for j = 1:nIter
    % Find where the maximum occurs on the current grid:
    [~ , ind] = max(abs(reshape(A, numel(A), 1)));
    % Indices of the maximum location
    [row , col] = ind2sub(size(A), ind);
    % Store the pivot location
    P(j,:) = [xpts(col), xpts(row)];
    % Update the function on the grid that we're approximating:
    A = A - A(:, col)*A(row, :)./A(row, col);
    % Update the function on the finer grid as well (for plotting):
    [~, ind] = max(abs(reshape(B, numel(B), 1)));
    [row, col] = ind2sub(size(B), ind);
    B = B - B(:, col)*B(row, :)./B(row, col);
    % Store the residual at the finer grid at each step:
    ff{j} = B;
    % Store the error after the current iteration:
    e(j) = norm(A);
end

end
