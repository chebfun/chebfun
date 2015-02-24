function movie(f,varargin)
%MOVIE play a chebfun2 in movie frames.
%
% movie(f) plays how the chebfun2 was constructed.
%
% This function is not a direct analogue of the movie function in MATLAB.
%
% Example 1:
% f = chebfun2(@(x,y) franke(x,y)); movie(f)
%
% Example 2:
% f = chebfun2(@(x,y) exp(-(x.^2+y.^2)/2)); movie(f)
%
% Example 3:
% f = chebfun2(@(x,y) cos(x.*y)); movie(f)

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.


% movie about the chebfun2 constructor.
home, close all
% Select the mode to play your chebfun2 in.
if nargin > 1
    if strcmpi(varargin{1},'speed') && nargin > 2
        mode = 'speed';
        speed=varargin{2};
    else
        mode = varargin{1};
    end
else
    mode = 'vslow';
end

if strcmpi(mode,'slow')
    ppan = 0.1; pmove=0.05; pfade =0.05; pbreak=.2; pwait=.5;psectionbreak=2;
elseif strcmpi(mode,'control')
    ppan = 0.1; pmove=0.05; pfade =0.05; pbreak=.2; pwait=.5;psectionbreak=inf;
elseif strcmpi(mode,'vslow')
    ppan = 0.2; pmove=0.1; pfade =0.1; pbreak=.4; pwait=1;psectionbreak=2;
elseif strcmpi(mode,'speed')
    c = 1/speed;
    ppan = 0.1*c; pmove=0.05*c; pfade =0.05*c; pbreak=.2*c; pwait=.5*c;psectionbreak=2*c;
else
    ppan = 0; pmove=0; pfade =0; pbreak=0; pwait=0;psectionbreak=0;
end

co = [0 0 1;0 1 0;1 0 0;0 1 1;1 0 1;1 1 0]; % color array.
pref = chebfunpref();
minsample = pref.minSamples; nextsample = 2^(ceil(log(minsample))+1)+1;

% Starting text
fprintf('This movie explains how a simplified version of chebfun2\nconstructor works on your function.\n \nHere goes...\n\n'),
txt = scribble('Your chebfun2 movie'); plot(txt), axis([-1 1 -1 1]), axis off
mypause(psectionbreak)
fprintf('Here is your chebfun2:\n\n'), clf

% Surface plot
% f = chebfun2(@(x,y) franke(x,y));
set(0,'DefaultAxesFontSize',14,'DefaultLineLineWidth',3)
sf=plot(f); axis equal, axis off,

% full pan around function.
sf=fullpan(sf,ppan); sf = tilt(sf,pmove,1); campos([0 0 10]), hold on
% slowing lower opaque
fade(1,.5,sf,pfade);mypause(psectionbreak)



% % % Start explaining the constructor % % %
fprintf(sprintf('We first sample on a %u x %u Chebyshev grid obtaining a matrix A.\n\n',minsample,minsample))
mypause(psectionbreak)

[xx yy]=meshgrid(chebpts(200));
B = feval(f,xx,yy);
% We first sample on a 9 by 9 grid.
[xx yy]=meshgrid(chebpts(minsample));
A = feval(f,xx,yy); scl = max(max(A));
zz = (max(max(A))-1e-1)*ones(minsample^2,1);
pts = plot3(xx(:),yy(:),zz,'Marker','.','MarkerSize',20,'Color','k','LineStyle','none');

% Take maximum on array, calculate residual and then take next maximum.
xpts = chebpts(minsample);
% Do first two steps of ACA.
P=zeros(2,2); g = @(x,y) f.feval(x,y);
for j = 1:min(2,length(f))
    [ infnorm , ind ]=max( abs ( reshape(A,numel(A),1) ) );
    [ row , col ]=ind2sub( size(A) , ind);
    P(j,:) = [xpts(col) xpts(row)];
    A = A - A( : , col )*A( row , : )./A(row,col);
    [ infnorm , ind ]=max( abs ( reshape(B,numel(B),1) ) );
    [ row , col ]=ind2sub( size(B) , ind);
    B = B - B( : , col )*B( row , : )./B(row,col);
    ff{j}=B;
    e(j) = norm(A);
end

fprintf('We then take entry with the largest absolute value.\n\n'),
mypause(psectionbreak)

colortext = {'blue','green'}; order = {'first','second'};
markers=[];
comet=newplot; t = 0:.01:1;  comets=[];
for jj= 1:min(2,length(f))
    val = feval(f,P(jj,1),P(jj,2)); %if (val > max(max(A))-1e-1), val = max(max(A))-1e-2; end % hack!
    mark = plot3(P(jj,1),P(jj,2),val,'Marker','.','Color',co(jj,:),'MarkerSize',40);
    if jj == 1
        fprintf('That''s the big blue dot and the first pivot.\n\n'),mypause(psectionbreak),
    else
        fprintf('That''s the green dot.\n\n');
    end
    mypause(psectionbreak)
    if jj==1
        fprintf('We take the entry''s value a, column u and row v...\n\n'), mypause(psectionbreak)
    end
    myline = mycomet(t,P(jj,:),comet,co(jj,:),pmove/5);
    mypause(psectionbreak)
    comets = [myline comets]; markers = [mark markers];
    
    if jj==1
        fprintf('and calculate the residual matrix A = A - uv^T/a.\n\n'), mypause(psectionbreak),
    end
    if jj==1
        fprintf(sprintf('The norm of the residual is %1.3e.\n\n',e(jj)));
    elseif jj==2
        fprintf(sprintf('After the second step the norm of the residual is %1.3e.\n\n',e(jj)))
    end
    
    if length(f)==1
        fprintf('The first step is complete.\n\n'), mypause(psectionbreak),
    elseif length(f) == 2 && jj==2
        fprintf('The first step is complete.\n\n'),mypause(psectionbreak),
    elseif length(f)>2 || (length(f) == 2 && jj==1)
        %         fprintf('The norm of the residual it was large.\n\n'),mypause(psectionbreak),
        if jj==1
            fprintf('We repeat for one more step on the residual matrix.\nAs before taking the largest absolute value.\n\n'),
            mypause(psectionbreak)
        end
    end
    delete(sf); rect = f.domain;
    sf = plot(chebfun2(ff{jj},rect)); alpha(sf,.5); colormap(jet)
    axis([rect,min(min(ff{jj}))-.1,scl+.1]);mypause(psectionbreak),
end
% delete(comets{:}), 
delete(markers), delete(pts);


if length(f) == 1
    % this was good enough, say so.
    fprintf('Your function is of numerical rank one.\n\n')
end

% Cannot take any more points so sample some more. Notice that the points
% interlace.
if length(f) > 2
    fprintf('Your function requires more sampling points to resolve.\n')
    mypause(psectionbreak)
    fprintf(sprintf('So we next sample it on a %u x %u Chebyshev grid.\n\n',nextsample,nextsample))
    mypause(psectionbreak)
    
    [xx yy]=meshgrid(chebpts(nextsample));
    A = feval(f,xx,yy);  scl = max(max(A)); sclinf = min(min(A));
    zz = (max(max(A))-1e-1)*ones(nextsample^2,1);
    pts = plot3(xx(:),yy(:),zz,'Marker','.','MarkerSize',20,'Color','k','LineStyle','none');
    
    [xx yy]=meshgrid(chebpts(200));
    B = feval(f,xx,yy);
    
    delete(sf), sf=plot(f); alpha(sf,.5); axis equal, axis off,
    
    if minsample ==9
        fprintf('The grids are nested for a little more efficiency.\nThe previous grid is shown in red.\n\n');
    end
    
    %         % fade from black to red to black.
    %         ptss=[];[xx yy]=meshgrid(chebpts(minsample)); zz = max(max(A))*ones(minsample^2,1);
    %         for t = -1:.1:-.01:.001:0:.001:.01:.1:1
    %             pts3 = plot3(xx(:),yy(:),zz,'Marker','.','MarkerSize',20,'Color',[1-abs(t),0,0],'LineStyle','none');
    %             ptss = [pts3 ptss];
    %             pause(pmove)
    %         end
    mypause(psectionbreak)
    fprintf('We repeat the same process just on a larger matrix...\n\n');
    % Same process again and now on a denser grid.
    xpts = chebpts(nextsample); [xx yy]=meshgrid(xpts);
    P=zeros(6,2);
    for j = 1:min(6,length(f))
        [ infnorm , ind ]=max( abs ( reshape(A,numel(A),1) ) );
        [ row , col ]=ind2sub( size(A) , ind);
        P(j,:) = [xpts(col) xpts(row)];
        A = A - A( : , col )*A( row , : )./A(row,col);
        [ infnorm , ind ]=max( abs ( reshape(B,numel(B),1) ) );
        [ row , col ]=ind2sub( size(B) , ind);
        B = B - B( : , col )*B( row , : )./B(row,col);
        ff{j}=B;
    end
    
    markers=[];
    comet=newplot; t = 0:.01:1;  comets=[];
    for jj= 1:min(6,length(f))
        val = feval(f,P(jj,1),P(jj,2)) + eps;
        mark = plot3(P(jj,1),P(jj,2),val,'Marker','.','Color',co(jj,:),'MarkerSize',40); hold on
        mypause(pbreak)
        myline = mycomet(t,P(jj,:),comet,co(jj,:),pmove/7);
        comets = [myline comets]; markers = [mark markers];
        delete(sf); rect = f.domain;
        sf = plot(chebfun2(ff{jj},rect));alpha(sf,.5); colormap(jet)
        axis([rect,min(min(ff{jj}))-.2,scl+.1]);mypause(psectionbreak),
    end
%         delete(comets{:}),
    delete(pts),
    %         delete(ptss),
    delete(sf)  % Remove everything so we can pan around.
    sf=plot(f);
    axis([rect,sclinf-.1,scl+.1]);
    axis equal, axis off,
else
    delete(sf)  % Remove everything so we can pan around.
    sf=plot(f);
    axis([rect,-scl-.1,scl+.1]);
    axis equal, axis off,
end
mypause(psectionbreak);

% If the function is very bad... need to tell people we take the first six
% slices only...
if length(f) > 6
    fprintf( sprintf('This process continues and your function was eventually sampled on a %u x %u Chebyshev grid. TODO:FIXME\n\n',-20))%length(f.fun2.C)));
end
if length(f)>1
    fprintf('The first stage is done.\n\n')
end
% Waterfall plot is what we have drawn.
water = waterfall(f,'-','nslices',min(length(f),6));

% tilt, pan, tilt
sf=tilt(sf,pmove,-1); sf=fullpan(sf,ppan); mypause(pbreak)
sf = tilt(sf,pmove,1); campos([0 0 10]), hold on

fprintf('To resolve your function we just need to resolve the selected columns and rows.\n\n'),mypause(2*pbreak);
fprintf('Each column and row is a function of one variable.\n\n'),mypause(psectionbreak);

delete(sf); mypause(psectionbreak)
% Skeleton version of ACA.
markers=[]; trun = min(6,length(f)); P = P(1:trun,:);
for j = 1:min(6,length(f))
    [xx yy]=meshgrid(P(:,1),chebpts(33));
    mark1 = plot3(xx,yy,scl+0*xx,'.k','MarkerSize',20);
    [xx yy]=meshgrid(chebpts(33),P(:,2));
    mark2 = plot3(xx,yy,scl+0*xx,'.k','MarkerSize',20);
    mypause(psectionbreak)
end
%     delete(mark1), delete(mark2),

fprintf('We use Chebfun to approximate it!\n\n'),mypause(psectionbreak);
fprintf('To check the columns and rows of your chebfun2 are resolved we can look\nat their Chebyshev coefficients.\n\n');


% fade surf plot away and delete.
%     fade(.75,0,sf,pfade);
delete(water), clf
%     if length(f) > 2, delete(markers), end
% Now draw chebpolyplot on the slices. Here are the columns.
cffs=newplot;
%     C = f.fun2.C; Ccfs = zeros(length(C),size(C,2)); col=[];
%     for jj = 1:min(6,length(f))
%         if pref2.mode
%             cc = chebpoly(C(:,jj));
%         else
%             cc = chebfft(C(:,jj));
%         end
%         cc(abs(cc)<eps)=eps;
%         Ccfs(:,jj) = cc;
%     end
%     Ccfs = abs(flipud(Ccfs)); %Ccfs = Ccfs(min(Ccfs.').'>0,:);
Ccfs = get(f.rows, 'coeffs');
Ccfs = log10(abs(flipud(Ccfs))); Ccfs(Ccfs==-inf)=NaN;
PivPos=f.pivotLocations;
cols=[];
for jj = 1:min(6,length(f))
    xx= PivPos(jj,1)*ones(length(Ccfs),1);
    yy = linspace(-1,1,length(Ccfs));
    col=plot3(xx,yy,Ccfs(:,jj),'Color',co(jj,:)); hold on,campos([0 0 10]),axis equal, axis off
    cols=[col cols];
end

% Now we do the rows.
%     R = f.fun2.R.'; Rcfs = zeros(length(R),size(R,2));
%     for jj = 1 : min(6,length(f))
%         if pref2.mode
%             rr = chebpoly(R(:,jj));
%         else
%             rr = chebfft(R(:,jj));
%         end
%         rr(abs(rr)<eps)=eps;
%         Rcfs(:,jj) = rr;
%     end
%     Rcfs = abs(flipud(Rcfs)); %Rcfs = Rcfs(min(Rcfs.').'>0,:);
Rcfs = get(f.cols, 'coeffs');
Rcfs = log10(abs(Rcfs));Rcfs(Rcfs==-inf)=NaN;
%     PivPos=f.fun2.PivPos;
rows=[];
for jj = 1:min(6,length(f))
    xx=linspace(-1,1,length(Rcfs));
    yy=PivPos(jj,1)*ones(length(Rcfs),1);
    row = plot3(xx,yy,Rcfs(:,jj),'Color',co(jj,:));
    rows=[row rows];
end
zlim([-18 0]), axis square;   % Fix axis zlimit before panning.

% tilt and then do some camorbits to get in the right position to see the
% coefficient decay.
sf=tilt(sf,pmove,-1); mypause(pbreak);

k=36;
for i=1:k
    camorbit(1,0,'data',[0 0 1]), camorbit(0,1,'data',[0 1 0]),mypause(pmove)
    drawnow
end
fprintf('This shows that the rows are resolved.\n\n')
mypause(psectionbreak)
for i=1:k
    camorbit(-1,0,'data',[0 0 1]), camorbit(0,-1,'data',[0 1 0]),mypause(pmove)
    drawnow
end
mypause(pbreak)
for i=1:k
    camorbit(-1.5,0,'data',[0 0 1]), camorbit(0,1,'data',[0 1 0]),mypause(pmove)
    drawnow
end
fprintf('This shows that the columns are resolved.\n\n')
mypause(psectionbreak)
for i=1:k
    camorbit(1.5,0,'data',[0 0 1]), camorbit(0,-1,'data',[0 1 0]),mypause(pmove)
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
%     delete(markers), delete(rows), delete(cols);
clf

fprintf('We store it away so you can play...\n\n'), mypause(pbreak)
% % % How are chebfun2 stored? % % %
txt = scribble('How is it stored?'); tt=plot(txt); axis([-1 1 -1 1]), axis off
mypause(psectionbreak), clf

% Make a plot that looks a bit like plot(f,'.-')
% f = chebfun2(@(x,y) franke(x,y));
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
fprintf('It''s stored as two chebfun quasi-matrices just like this:\n\n'), mypause(pbreak)
for t = 0:.1:1
    delete(cols),delete(rows), delete(markers)
    cols=[];
    for jj= 1:min(6,length(f))
        st = t*(cc(jj)-crosses(jj,1));
        col = line([crosses(jj,1)+st crosses(jj,1)+st].',[-1 1],'Color',co(jj,:)); hold on,
        cols = [col cols];
    end
    rows=[];
    for jj= 1:min(6,length(f))
        st = t*(rr(jj)-crosses(jj,2));
        row = line([-1+1.5*t 1+1.5*t],[crosses(jj,2)+st crosses(jj,2)+st].','Color',co(jj,:));
        rows = [row rows];
    end
    markers=[];
    for jj= 1:min(6,length(f))
        stx = t*(dd(jj)-P(jj,1)); sty = t*(rr(jj)-P(jj,2));
        mark = plot(P(jj,1)+stx,P(jj,2)+sty,'Marker','.','MarkerSize',40,'Color',co(jj,:));
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


function h=fullpan(h,p)
k=36; %unit = pi/k;
for i=1:k
    camorbit(10,0,'data',[0 0 1])
    mypause(p)
    drawnow
end
end

function h=tilt(h,p,whichway)
k=36;
s = whichway;
for i=1:k
    camorbit(0,s*1.665,'data',[0 1 1])
    camorbit(s*1.04,0,'data',[0 0 1])
    mypause(p)
    drawnow
end
end

function fade(alpha1,alpha2,h,p)
for a = alpha1:-.1:alpha2
    alpha(h,a); mypause(p);
end
end

function h = mycomet(t,pivot,h,color,p)
x=zeros(1,length(t));
xeven = pivot(1) + t(1:2:end)*(1-pivot(1));
xodd = pivot(1) + t(2:2:end)*(-1-pivot(1));
x(1:2:end)=xeven; x(2:2:end)=xodd;
y = pivot(2)*ones(length(t),1);

% if verLessThan('matlab', '8.4')
body = line('parent',h,'color',color,'linestyle','-','LineWidth',3,'erase','none', ...
    'xdata',[],'ydata',[]);
% else
%     body = animatedline('parent',h,'color',color,'linestyle','-','LineWidth',3,'erase','none', ...
%         'xdata',[],'ydata',[]);
% end
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

h={[body body2]};

end


function mypause(n)
if n==inf
    pause;
else
    pause(n)
end
end
