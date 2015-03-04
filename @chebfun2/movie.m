function movie(f, varargin)
%MOVIE play a a movie about the CHEBFUN2 constructor
%
% MOVIE(F), where F is a CHEBFUN2, plays a movie about how the CHEBFUN2
% constructor constructed F.
%
% MOVIE(F, MODE) speeds up or slow the movie, depending on the value of MODE.
% The options for the value of MODE are as follows:
%
%   MODE = S, where S is a numerical value, plays the movie S times faster than
%   the regular speed.
%   
%   MODE = 'SLOW': TODO: Describe
%
%   MODE = 'VSLOW': TODO: Describe
%
% Examples:
%   f = chebfun2(@(x,y) franke(x,y)); movie(f)
%   f = chebfun2(@(x,y) exp(-(x.^2+y.^2)/2)); movie(f, 2)
%   f = chebfun2(@(x,y) cos(x.*y)); movie(f, 'slow')
%
% Note: This function is not a direct analogue of the movie function in MATLAB.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Tidy up before we start the movie:
% home, close all

%% Setup the options for playing the movie, in particular, the speed of it.

% Select the mode to play your chebfun2 in.
if nargin > 1
    if strcmpi(varargin{1}, 'speed') && nargin > 2
        mode = 'speed';
        speed=varargin{2};
    else
        mode = varargin{1};
    end
else
    mode = 'vslow';
end

% Setup some parameters for the screening of the movie. The parameters govern
% the following attributes:
%   PPAN:
%   PMOVE:
%   PFADE: Pause between steps while fading
%   PBREAK:
%   PSECTIONBREAK:
if strcmpi(mode, 'slow')
    ppan   = 0.1;
    pmove  = 0.05;
    pfade  = 0.05;
    pbreak = 0.2;
    psectionbreak = 2;
elseif strcmpi(mode, 'control')
    ppan   = 0.1;
    pmove  = 0.05;
    pfade  = 0.05;
    pbreak = 0.2;
    psectionbreak = inf;
elseif strcmpi(mode, 'vslow')
    ppan   = 0.2;
    pmove  = 0.1;
    pfade  = 0.1;
    pbreak = 0.4;
    psectionbreak = 2;
elseif strcmpi(mode, 'speed')
    c = 1/speed;
    ppan   = 0.1*c;
    pmove  = 0.05*c;
    pfade  = 0.05*c;
    pbreak = 0.2*c;
    psectionbreak = 2*c;
else
    ppan   = 0;
    pmove  = 0;
    pfade  = 0;
    pbreak = 0;
    psectionbreak = 0;
end

% Colors used for plotting lines on the plot. This is useful for old versions of
% Matlab compatability.
co = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];

% For the movie, we only want to look at 9 and 17 point grids (although the
% current preference of MINSAMPLES might be different), as otherwise, it becomes
% too cluttered.
pref = chebfunpref();
minsample = 9;pref.minSamples;
nextsample = 17;2^(ceil(log(minsample)) + 2) + 1;

%% Print starting text
fprintf('This movie explains how a simplified version of the chebfun2\n')
fprintf('constructor works on your function.\n\n')
fprintf('Here it goes...\n\n')
% Scribble to the plot
txt = scribble('Your chebfun2 movie');
plot(txt)
axis([-1 1 -1 1]), axis off
mypause(psectionbreak)

%% Surface plot of the input CHEBFUN2, followed by panning around the function
fprintf('This is your chebfun2:\n\n'), clf
set(0, 'DefaultAxesFontSize', 14, 'DefaultLineLineWidth', 3)
sf = plot(f);
axis equal, axis off

% Full pan around function:
sf = fullpan(sf, ppan);
% Tilt to a view from above:
sf = tilt(sf, pmove, 1);
% Slowly go from 100% to 50% opaque:
fade(1, .5, sf, pfade);
% Pause before next step:
mypause(psectionbreak)

%% Start explaining the constructor
fprintf(['We first sample on a %u x %u Chebyshev grid obtaining ' ...
    'a matrix A.\n\n'], minsample, minsample)
mypause(psectionbreak)

% Sample the function on a finer Chebyshev grid so that we can compute (and
% plot) the residual after taking one step.
[xx, yy] = meshgrid(chebpts(200));
B = feval(f, xx, yy);

% We first sample on a 17x17 grid.
[xx, yy] = meshgrid(chebpts(minsample));
A = feval(f, xx, yy); scl = max(max(A));
zz = (max(max(A)) - 1e-1)*ones(minsample^2, 1);
hold on
pts = plot3(xx(:), yy(:), zz, 'Marker', '.', 'MarkerSize', 20, ...
    'Color', 'k', 'LineStyle', 'none');

% Take maximum on array, calculate residual and then take next maximum.
xpts = chebpts(minsample);
% Store the pivot locations
P = zeros(2, 2);
for j = 1:min(2, length(f))
    % Find where the maximum occurs on the current grid:
    [infnorm , ind] = max(abs(reshape(A, numel(A), 1)));
    % Indices of the maximum location
    [ row , col ] = ind2sub(size(A) , ind);
    % Store the pivot location
    P(j,:) = [xpts(col), xpts(row)];
    % Update the function on the grid that we're approximating:
    A = A - A(:, col)*A(row, :)./A(row, col);
    % Update the function on the finer grid as well (for plotting):
    [infnorm, ind]=max(abs(reshape(B, numel(B), 1)));
    [row, col]=ind2sub(size(B), ind);
    B = B - B(:, col )*B(row, : )./B(row,col);
    % Store the residual at the finer grid at each step:
    ff{j}=B;
    % Store the error after the current iteration:
    e(j) = norm(A);
end

fprintf('We then take entry with the largest absolute value.\n\n'),
mypause(psectionbreak)

markers=[];
comet=newplot; t = 0:.01:1;  comets=[];
for jj= 1:min(2, length(f))
    val = feval(f, P(jj,1), P(jj,2));
    mark = plot3(P(jj,1), P(jj,2), val, 'Marker', '.', ...
        'Color', co(jj,:), 'MarkerSize', 40);
    if ( jj == 1 )
        fprintf('That''s the big blue dot and the first pivot.\n\n')
        mypause(psectionbreak)
    else
        fprintf('That''s the red dot.\n\n')
    end
    mypause(psectionbreak)
    if ( jj==1 )
        fprintf('We take the entry''s value a, column u and row v...\n\n')
        mypause(psectionbreak)
    end
    myline = mycomet(t, P(jj,:), comet,co(jj,:), pmove/5);
    mypause(psectionbreak)
    comets = [myline comets];
    markers = [mark markers];
    
    if ( jj==1 )
        fprintf('and calculate the residual matrix A = A - uv^T/a.\n\n'), mypause(psectionbreak),
    end
    if ( jj==1 )
        fprintf(sprintf('The norm of the residual is %1.3e.\n\n',e(jj)));
    elseif ( jj==2 )
        fprintf(sprintf('After the second step the norm of the residual is %1.3e.\n\n',e(jj)))
    end
    
    if ( length(f) == 1 )
        fprintf('The first step is complete.\n\n')
        mypause(psectionbreak),
    elseif ( length(f) == 2 && jj == 2 )
        fprintf('The first step is complete.\n\n')
        mypause(psectionbreak),
    elseif ( length(f) > 2 || (length(f) == 2 && jj ==1 ) )
        if ( jj == 1 )
            fprintf('We repeat for one more step on the residual matrix.\n')
            fprintf('As before taking the largest absolute value.\n\n')
            mypause(psectionbreak)
        end
    end
    delete(sf); rect = f.domain;
    sf = plot(chebfun2(ff{jj}, rect));
    alpha(sf,.5);
    axis([rect,min(min(ff{jj})) - .1, scl + .1]);
    mypause(psectionbreak),
end
% TODO: Fails for rank1 chebfun2
delete(comets{1})
if length(f) > 1
    delete(comets{2})
end
delete(markers), delete(pts);


if ( length(f) == 1 )
    % this was good enough, say so.
    fprintf('Your function is of numerical rank one.\n\n')
end

%% We did not resolve the function with two pivot points, so take more.
if ( length(f) > 2 )
    fprintf('Your function requires more sampling points to resolve.\n')
    mypause(psectionbreak)
    fprintf(sprintf('So we next sample it on a %u x %u Chebyshev grid.\n\n', ...
        nextsample, nextsample))
    mypause(psectionbreak)
 
    [xx, yy] = meshgrid(chebpts(nextsample));
    A = feval(f, xx, yy);  scl = max(max(A)); sclinf = min(min(A));
    zz = (max(max(A)) - 1e-1)*ones(nextsample^2, 1);
    pts = plot3(xx(:), yy(:), zz, 'Marker', '.', 'MarkerSize', 20,...
       'Color','k','LineStyle','none');
 
    [xx, yy]=meshgrid(chebpts(200));
    B = feval(f, xx, yy);
 
    delete(sf), sf=plot(f); alpha(sf,.5); axis equal, axis off
    
    mypause(psectionbreak)
    fprintf('We repeat the same process as before, just on a larger matrix...\n\n');
    % Same process again and now on a denser grid.
    xpts = chebpts(nextsample);
    [xx, yy] = meshgrid(xpts);
    P = zeros(6, 2);
    for j = 1:min(6, length(f))
        [infnorm , ind] = max( abs ( reshape(A, numel(A), 1) ) );
        [ row , col ] = ind2sub( size(A) , ind);
        P(j,:) = [xpts(col) xpts(row)];
        A = A - A( : , col )*A( row , : )./A(row,col);
        [infnorm , ind] = max( abs ( reshape(B,numel(B),1) ) );
        [row , col] = ind2sub(size(B), ind);
        B = B - B( : , col )*B( row , : )./B(row,col);
        ff{j}=B;
    end
    
    markers=[];
    comet=newplot; t = 0:.01:1;  comets=[];
    for jj= 1:min(6, length(f))
        val = feval(f,P(jj,1),P(jj,2)) + eps;
        mark = plot3(P(jj,1), P(jj,2), val, 'Marker', '.', ...
            'Color', co(jj,:), 'MarkerSize', 40);
        hold on
        mypause(pbreak)
        myline = mycomet(t, P(jj,:), comet, co(jj,:), pmove/7);
        comets = [myline comets];
        markers = [mark markers];
        delete(sf);
        rect = f.domain;
        sf = plot(chebfun2(ff{jj}, rect)); alpha(sf,.5);
        axis([rect,min(min(ff{jj})) - .2, scl + .1]);
        mypause(psectionbreak),
    end
    for comCounter=1:length(comets)
        delete(comets{comCounter})
    end
    delete(pts)
    delete(sf)  % Remove everything so we can pan around.
    sf=plot(f);
    axis([rect, sclinf - .1, scl + .1]);
    axis equal, axis off
else
    delete(sf)  % Remove everything so we can pan around.
    sf=plot(f);
    axis([rect, -scl - .1, scl + .1]);
    axis equal, axis off
end
mypause(psectionbreak);

% If the function is very bad... need to tell people we take the first six
% slices only...
if ( length(f) > 6 )
    fprintf('This process continues and your function was eventually sampled on a %u x %u Chebyshev grid.\n\n', length(f.rows(:,1)));
elseif ( length(f) > 1 )
    fprintf('The first stage is done.\n\n')
end

%% Waterfall plot is what we have drawn.
water = waterfall(f, '-', 'nslices', min(length(f), 6));

% tilt, pan, tilt
sf = tilt(sf, pmove, -1); sf = fullpan(sf, ppan); mypause(pbreak)
sf = tilt(sf, pmove, 1); campos([0 0 10]), hold on

fprintf('To resolve your function we just need to resolve the selected columns and rows.\n\n')
mypause(2*pbreak);
fprintf('Each column and row is a function of one variable.\n\n')
mypause(psectionbreak);

delete(sf);
mypause(psectionbreak)

%% Skeleton version of ACA.
markers=[]; trun = min(6,length(f)); P = P(1:trun,:);
for j = 1:min(6, length(f))
    [xx, yy] = meshgrid(P(:,1), chebpts(33));
    mark1 = plot3(xx,yy,scl+0*xx,'.k','MarkerSize',20);
    [xx, yy] = meshgrid(chebpts(33), P(:,2));
    mark2 = plot3(xx, yy, scl+0*xx, '.k', 'MarkerSize', 20);
    mypause(psectionbreak)
end
%     delete(mark1), delete(mark2),
delete(water), clf

%% Do chebpolyplots of all the columns and rows
fprintf('We use Chebfun to approximate it!\n\n'),mypause(psectionbreak);
fprintf('To check the columns and rows of your chebfun2 are resolved we can look\nat their Chebyshev coefficients.\n\n');

% Draw chebpolyplot on the slices. Here are the columns.
Ccfs = get(f.cols, 'coeffs');
% Only want the first six columns
if length(f) > 6
    Ccfs(:, 7:length(f)) = [];
end

Ccfs = log10(abs(flipud(Ccfs))); Ccfs(Ccfs==-inf)=log10(eps);
PivPos = f.pivotLocations;
cols = [];
for jj = 1:min(6,length(f))
    xx = PivPos(jj,1)*ones(length(Ccfs), 1);
    yy = linspace(-1, 1, length(Ccfs));
    col=plot3(xx, yy, Ccfs(:,jj), 'Color', co(jj,:));
    hold on
    campos([0 0 10])
    axis equal, axis off
    cols=[col cols];
end

% Now we do the rows.
Rcfs = get(f.rows, 'coeffs');
% Only want the first six rows:
if length(f) > 6
    Rcfs(:, 7:length(f)) = [];
end
Rcfs = log10(abs(Rcfs));Rcfs(Rcfs==-inf)=log10(eps);

rows=[];
for jj = 1:min(6,length(f))
    xx=linspace(-1,1,length(Rcfs));
    yy=PivPos(jj,1)*ones(length(Rcfs),1);
    row = plot3(xx,yy,Rcfs(:,jj),'Color',co(jj,:));
    rows=[row rows];
end

%% Tilt and camorbit to show the coefficients decaying
% Fix axis zlimit before panning.
zlim([min(min(min(Ccfs)), min(min(Rcfs))) 0]), axis square;
shg
sf = tilt(sf, pmove, -1);
mypause(pbreak);

k = 36;
for i = 1:k
    camorbit(1, 0, 'data', [0 0 1])
    camorbit(0, 1 ,'data', [0 1 0])
    mypause(pmove)
    drawnow
end
fprintf('This shows that the rows are resolved.\n\n')
mypause(psectionbreak)
for i = 1:k
    camorbit(-1, 0, 'data', [0 0 1])
    camorbit(0, -1, 'data', [0 1 0])
    mypause(pmove)
    drawnow
end
mypause(pbreak)
for i = 1:k
    camorbit(-1.5, 0, 'data', [0 0 1])
    camorbit(0, 1, 'data', [0 1 0]), mypause(pmove)
    drawnow
end
fprintf('This shows that the columns are resolved.\n\n')
mypause(psectionbreak)
for i = 1:k
    camorbit(1.5, 0, 'data', [0 0 1])
    camorbit(0, -1, 'data', [0 1 0])
    mypause(pmove)
    drawnow
end
mypause(pbreak), sf=tilt(sf,pmove,1);


% Put back the pivot locations.
P = PivPos; U = f.pivotValues; markers=[];
for jj= 1:min(6,length(f))
    val = feval(f,P(jj,1),P(jj,2)) + eps;
    mark = plot3(P(jj,1),P(jj,2),val,'Color',co(jj,:),'Marker','.','MarkerSize',40);
    markers = [mark markers];
end
campos([0 0 10]);

mypause(pbreak),
clf

%% Show how the chebfun2 is stored
fprintf('We store it away so you can play...\n\n'), mypause(pbreak)
% % % How are chebfun2 stored? % % %
txt = scribble('How is it stored?');
plot(txt); axis([-1 1 -1 1]), axis off
mypause(psectionbreak), clf

% Make a plot that looks a bit like plot(f, '.-')
crosses = PivPos;
cols = line([crosses(:,1) crosses(:,1)].',[-1 1]); hold on,
rows = line([-1 1],[crosses(:,2),crosses(:,2)].');
mx = 2.3; axis(mx*[-1 1 -1 1]), axis square, axis off

% Again plot pivot positions.
P = PivPos; U = f.pivotValues; markers=[];
for jj= 1:min(6,length(f))
    mark = plot(P(jj,1),P(jj,2),'Marker','.','MarkerSize',40,'Color',co(jj,:));
    markers = [mark markers];
end

% Slowly move columns and rows to low rank storage structure.
if length(f) > 1
    cc = linspace(-1.75,-1,min(length(f),6));
    rr = linspace(.25,1,min(length(f),6)); rr=rr(end:-1:1);
    dd = linspace(-.6,.15,min(length(f),6));
else
    cc = -1; rr = 1; dd=-.6;
end
fprintf('It''s stored as two chebfun quasi-matrices just like this:\n\n')
mypause(pbreak)
for t = 0:.1:1
    delete(cols),delete(rows), delete(markers)
    cols=[];
    for jj= 1:min(6,length(f))
        st = t*(cc(jj)-crosses(jj,1));
        col = line([crosses(jj,1)+st crosses(jj,1)+st].', [-1 1], ...
            'Color',co(jj,:));
        hold on
        cols = [col cols];
    end
    rows=[];
    for jj= 1:min(6,length(f))
        st = t*(rr(jj)-crosses(jj,2));
        row = line([-1+1.5*t 1+1.5*t],[crosses(jj,2)+st crosses(jj,2)+st].', ...
            'Color',co(jj,:));
        rows = [row rows];
    end
    markers=[];
    for jj= 1:min(6,length(f))
        stx = t*(dd(jj)-P(jj,1)); sty = t*(rr(jj)-P(jj,2));
        mark = plot(P(jj,1) + stx, P(jj,2) + sty, 'Marker', '.', ...
            'MarkerSize', 40, 'Color', co(jj,:));
        markers = [mark markers];
    end
    mypause(pbreak);
end
% Create textbox
annotation(gcf,'textbox',...
    [0.143441860465116 0.201923076923077 0.764988372093023 0.0913461538461542],...
    'String',{'f(x,y) = C(y)         U          R(x)'},...
    'FontSize',20,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'LineWidth',3);
end


function h = fullpan(h, p)
% Set the camera view angle to manual mode, so that the size of the plot
% doesn't keep changing as we rotate
set(gca, 'cameraviewanglemode', 'manual')
axis vis3d
k = 36;
for i = 1:k
    % Spin the camera around the z-axis (the third and fourth arguments to
    % camorbit specify what axis we want to spin around)
    camorbit(10, 0, 'data', [0 0 1])
    mypause(p)
    drawnow
end
end

function h = tilt(h, p, whichway)
set(gca,'cameraviewanglemode','auto')
k = 36;
s = whichway;
for i = 1:k
    camorbit(0, s*1.665, 'data', [0 1 1])
    camorbit(s*1.04, 0, 'data', [0 0 1])
    mypause(p)
    drawnow
end
end

function fade(alpha1, alpha2, h, p)
%FADE
for a = alpha1:-.1:alpha2
    alpha(h, a);
    mypause(p);
end
end

function h = mycomet(t, pivot, h, color, p)
x=zeros(1, length(t));
xeven = pivot(1) + t(1:2:end)*(1 - pivot(1));
xodd = pivot(1) + t(2:2:end)*(-1 - pivot(1));
x(1:2:end)=xeven;
x(2:2:end)=xodd;
y = pivot(2)*ones(length(t), 1);

if verLessThan('matlab', '8.4')
    body = line('parent',h,'color',color,'linestyle','-','LineWidth',3,'erase','none', ...
        'xdata',[],'ydata',[]);
    for jj = 1
        set(body,'xdata',x(jj),'ydata',y(jj)), drawnow, hold on
    end
    for jj = 2:length(x)
        j = jj-1:jj;
        set(body,'xdata',x(j),'ydata',y(j)), drawnow
        mypause(p)
    end
    
    body2 = line('parent',h,'color',color,'linestyle','-','LineWidth',2,'erase','none', ...
        'xdata',[],'ydata',[]);
    x = pivot(1)*ones(length(t),1);
    yeven = pivot(2) + t(1:2:end)*(1-pivot(2));
    yodd = pivot(2) + t(2:2:end)*(-1-pivot(2));
    y(1:2:end)=yeven; y(2:2:end)=yodd;
    for jj = 1
        set(body2,'xdata',x(jj),'ydata',y(jj)), drawnow
    end
    for jj = 2:length(x)
        j = jj-1:jj;
        set(body2,'xdata',x(j),'ydata',y(j)), drawnow
        mypause(p)
    end
else
%     body = animatedline('parent',h,'color',color,'linestyle','-','LineWidth',3,'erase','none', ...
%         'xdata',[],'ydata',[]);
    body = animatedline;
    body2 = animatedline;
    for k = 1:length(x);
        addpoints(body, x(k), y(k))
        drawnow update
    end
    
    x = pivot(1)*ones(length(t),1);
    yeven = pivot(2) + t(1:2:end)*(1-pivot(2));
    yodd = pivot(2) + t(2:2:end)*(-1-pivot(2));
    y(1:2:end)=yeven; y(2:2:end)=yodd;
    for k = 1:length(x);
        addpoints(body2, x(k), y(k))
        drawnow update
    end
    
end


h={[body body2]};

end


function mypause(n)
if ( n == inf ) 
    pause
else
    pause(n)
end
end
