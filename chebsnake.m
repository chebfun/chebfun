function chebsnake(nodes,alpha)
%CHEBSNAKE   Chebfun snake game.
%   CHEBSNAKE() plays a twist on a classic game where you must feed the snake
%   with more and more interpolation nodes, but avoid that it hits the boundary
%   or itself! Use the arrow keys to control the snake. Any other key will quit
%   the game.
%
%   CHEBSNAKE(NODES) allows one to change the interpolation type. The default
%   type 'cheby' is polynomial interpolation in Chebyshev points. Other types
%   are polynomial interpolation in equispaced points ('equi') and
%   Floater-Hormann rational interpolation in equispaced points ('fh'). The blue
%   dots on the snake indicate the interpolated function values.
%
%   CHEBSNAKE(NODES, ALPHA) allows to change the initial game speed by a factor
%   ALPHA > 0, ALPHA > 1 increases the game speed, ALPHA < 1 decreases it
%   (default = 1).
%
%   To prevent you from neglecting your actual work, the game speed increases
%   with the total number of achieved points...
%
% See also CHEBTUNE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get some constants right:
W = warning;
if ( nargin < 2 )
    alpha = 1;
end
if ( nargin > 0 && strcmp(nodes, 'equi') )
    nodes = 0;
    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
elseif ( nargin > 0 && strcmp(nodes, 'fh') ) 
    nodes = 2;
else
    nodes = 1;
end
LW = 'LineWidth'; lw = 2;
MFC= 'MarkerFaceColor';
res = 0.15; 
len = 5; 
dom = domain([-1, 1]); 
d = 1;
food = @() res*(round((1.8*rand-.9)/res)+1i*round((1.8*rand-.9)/res));
pause on

% keyboard interaction
figure('KeyPressFcn', @keypress);
    function keypress(~, evnt)
        dold = d;
        switch evnt.Key
            case 'leftarrow',   d = -1;  % left
            case 'rightarrow',  d = 1;   % right
            case 'downarrow',   d = -1i; % down
            case 'uparrow',     d = 1i;  % up
            otherwise,          d = 0;   % quit
        end;
        if ( d == -dold )
            d = dold;
        end
    end

lvl = 1;
pts = 0;
fails = 0;                          % fail counter (no food eaten)
failmax = 5;                        % number of consecutive fails before quit

lv = chebfun(@(x) exp(x), [0 pi]);
dd = chebfun(@(x) exp(-x), [-pi 0]);
grb = chebfun(@(x) 20*cos(x), [0 2*pi]);
kld = chebfun(@(x) exp(-.5*x).*(2+cos(10*x)), [-pi 0]);

while ( d ~= 0 ) % until quit
    d = 1;
    clf
    
    s = linspace(res*(1-len), 0, len) + 1i*eps;
    hs1 = plot(s(1:end-1), 'b-', LW, lw); hold on
    hs2 = [ plot(s(1:end-2), 'bo', MFC, 'b'), ...
            plot(s(end-1), 'bo', LW, lw, 'MarkerSize', 8) ];
    f = food();
    hf = plot(real(f), imag(f), 'ro', 'MarkerSize', 10);
    if (~rem(pts + 1, 30) && pts )
        set(hf, MFC, 'r', 'Color', [0, 0.6, 0]);
    else
        set(hf, MFC, [0, 0.6, 0], 'Color', 'r');
    end
    ht = plot(8, 0); axis square;       % dummy handle
    set(gca,'XTick', [], 'YTick', []);
    title('Control the snake with arrow keys. Quit with any other key.');
    axis([-1, 1, -1, 1]); 
    shg; 
    pause(0.3);
    
    % Ready, set, go!
    t = 1;                              % convex factor for nodes
    go = chebsnakePlotChebfun(.7*scribble('ready?'), 'r', LW, lw);
    shg; 
    pause(0.3);
    delete(go)
    
    go = chebsnakePlotChebfun(.4*scribble('go!'), 'r', LW, lw);
    shg; 
    pause(0.3); 
    delete(go)
    
    tic
    while ( d ~= 0 )                    % until game over or quit
        
        t = t + .2*alpha;
        if ( t > 1 )
            t = 0; 
            dr = res*d;
            s = [ s(2:end), s(end) + dr ];
            if ( length(s) < len + pts )
                s = [ s(2), s ];  %#ok<AGROW>
            end
        end
        
        y = (1 - t)*s(1:end-1) + t*s(2:end);
        if ( nodes == 1 )
            c = chebfun(y.');
        elseif ( nodes == 2 )
            fhd = min(ceil(0.4*sqrt(length(y))), 4);
            yy = linspace(-1, 1, 5*length(y)).';
            xx = linspace(-1, 1, length(y)).';
            ww = fhweights(length(y)-1, fhd);
            c = bary(yy, y.', xx, ww);
        elseif ( nodes == 0 )
            yy = linspace(-1, 1, length(y)).';
            c = polyfit(yy, y, length(y) - 1, dom);
        end
        for k = 1:numel(hs1)
            delete(hs1(k));
        end

        hs1 = chebsnakePlotChebfun(c, 'b-', LW, lw);
        delete(hs2);
        hs2 = [plot(y(1:end-1), 'bo', MFC, 'b'), ...
               plot(y(end), 'bo', LW, lw, 'MarkerSize', 8)];
        set(gca, 'xLim', [-1 1], 'yLim', [-1 1]);
        shg
        pause(max(0.01, 0.03 - toc) / alpha); 
        
        tic
        % check if the snake hits itself or the boundary
        if ( (max(abs([real(y(end)), imag(y(end))])) > 1) || ...
                (min(abs(y(end)-y(1:end-1))) < res/2) )
            ht = chebsnakePlotChebfun(.8*scribble('game over'), 'r', LW, lw);
            chebtune(dd, .5);
            shg
            pause(1)
            fails = fails + 1;
            if ( fails > failmax )
                d = 0; 
            end
            break
        end
        
        if ( abs(y(end) - f) < res/2 ) % snake eats food ?
            
            pts = pts + 1;
            chebtune(grb, .5);
            if ( ~rem(pts, 10) )
                lvl = lvl + 1;
                alpha = alpha * 1.1;
                chebtune(10*chebpoly(pts));
            end
            
            if ( ~rem(pts, 30) )
                chebtune(lv, 1);
                fails = fails - 1;
                up = chebsnakePlotChebfun(.8*scribble('1 up!'), 'r', LW, lw);
                shg
                pause(1)
                delete(up)
            end
            
            title(['Points : ' num2str(pts) '       Level : ' num2str(lvl) ...
                '       Lives: ' num2str(failmax - fails)], 'color', 'k');
            
            f = food();
            while ( any( abs(f - y) < res/2) )
                f = food();
            end
            set(hf, 'XData', real(f), 'YData', imag(f));
            if ( ~rem(pts + 1, 30) )
                set(hf, MFC, 'r', 'Color', [0,0.6,0]);
            else
                set(hf, MFC, [0, 0.6, 0], 'Color', 'r');
            end
            
        end
    end
    
    for k = 1:numel(ht)
        delete(ht(k));
    end
    
end

chebsnakePlotChebfun(.8*scribble('goodbye'), 'r', LW, lw); 
chebtune(kld, 1);
shg
pause(1)
close(gcf)
warning(W)

    function w = fhweights(n, fhd) 
    % Weights for Floater-Hormann interpolation
        w = zeros(1, n+1);
        for l = 0:n
            ji = max(l - fhd, 0);
            jf = min(l, n - fhd);
            sumcoeff = zeros(jf - ji + 1, 1);
            for i = ji:jf
                sumcoeff(i-ji+1) = nchoosek(fhd, l - i);
            end
            w(l+1) = (-1)^(l-fhd)*sum(sumcoeff);
        end
    end

end

function h = chebsnakePlotChebfun(varargin)
%CHEBSNAKEPLOT   Plot a CHEBFUN used in CHEBSNAKE.
%   This function is just a wrapper for CHEBFUN PLOT which gets all the plot
%   handles and assembles them into a numeric vector, making them easier to
%   free later.
    [hl, hp, hj, hd] = plot(varargin{:});
    h = [hl ; hp ; hj ; hd];
end
