function chebsnake2( f, nodes, alfa )
%CHEBSNAKE2   CHEBFUN2 snake game on a surface.
%   This is the CHEBFUN2 analogue of Chebfun's CHEBSNAKE.
%
%   CHEBSNAKE2 Feed the snake with more and more 1D interpolation nodes, but
%   avoid that it hits the boundary or itself! Use the arrow keys to control the
%   snake. Any other key will quit the game.
%
%   CHEBSNAKE2(F) allows one to specify the CHEBFUN2 object F on which snake
%   will live (default is chebfun2(@(x, y) 2-x.^2 - y.^2)).
%
%   CHEBSNAKE2(F, NODES) allows one to change the interpolation nodes and type.
%   The default type 'cheby' is polynomial interpolation in Chebyshev points.
%   Other types are polynomial interpolation in equispaced points ('equi') and
%   Floater-Hormann rational interpolation in equispaced points ('fh'). The blue
%   dots on the snake indicate the interpolated function values.
%
%   CHEBSNAKE2(F, NODES, ALFA) allows to change the initial game speed by a
%   factor alfa > 0, alfa > 1 increases the game speed, alfa < 1
%   decreases it (default = 1).
%
%   To prevent you from neglecting your actual work, the game speed
%   increases with the total number of achieved points.
%
% See also CHEBSNAKE, CHEBTUNE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: This code needs to be cleaned up. A lot.

% Get some constants right
W = warning;
if ( nargin < 3 )
    alfa = 1;
end
if ( nargin < 2 )
    nodes = 1;
end
if ( (nargin > 0 && ~isa(f, 'chebfun2')) || nargin == 0 || isempty(f) )
    f = chebfun2(@(x, y) 2-x.^2 - y.^2);
end
if ( nargin > 1 && strcmp(nodes, 'equi') )
    nodes = 0;  warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
elseif ( nargin > 1 && strcmp(nodes, 'fh'))
    nodes = 2;
else
    nodes = 1;
end

LW = 'LineWidth'; lw = 2;
MS = 'MarkerSize'; ms = 10;
res = 0.15; 
len = 5; 
dom = domain([-1, 1]); 
d = 1;
food = @() res*(round((1.8*rand-.9)/res)+1i*round((1.8*rand-.9)/res));
pause on

minf = min2( f );
if ( minf < 0 )
    f = f + abs( min2( f ) );
end
maxf = max2( f );
if maxf == 0
    maxf = 1;
end

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
fails = 0;
hold on
view(0, 4)
failmax = 5; % number of consecutive fails before quit
lv = chebfun(@(x) exp(x), [0 pi]);
dd = chebfun(@(x) exp(-x), [-pi 0]);
grb = chebfun(@(x) 20*cos(x), [0 2*pi]);
kld = chebfun(@(x) exp(-.5*x).*(2+cos(10*x)), [-pi 0]);

while (d ~= 0), % until quit
    d = 1;
    
    clf
    plot(f)
    colormap hsv, alpha(.2)
    hold on
    s = linspace(res*(1-len), 0, len) + 1i*eps;
    hs1 = plot3(real(s(1:end-1)), imag(s(1:end-1)), ...
        f(real(s(1:end-1)), imag(s(1:end-1))), 'b-', LW, lw); hold on
    hs2 = [plot3(real(s(1:end-2)), imag(s(1:end-2)), f(real(s(1:end-2)), ...
        imag(s(1:end-2))), 'bo', 'MarkerFaceColor', 'b', LW, lw) ...
        plot3(real(s(end-1)), imag(s(end-1)), f(real(s(end-1)), ...
        imag(s(end-1))), 'bo', LW, lw)];
    hs1s = plot(s(1:end-1), 'k-', LW, lw/2);
    hs2s = plot(s(1:end-1), 'k.', LW, lw/2);
    fd = food();
    hf = plot3(real(fd), imag(fd), f(real(fd), imag(fd)), 'ro', MS, ms);
    if ( ~rem(pts+1, 30) && pts )
        set(hf, 'MarkerFaceColor', 'r', 'Color', [0, 0.6, 0]);
    else
        set(hf, 'MarkerFaceColor', [0, 0.6, 0], 'Color', 'r');
    end
    hfs = plot(real(fd), imag(fd), 'ko', MS, ms/2, 'MarkerFaceColor', 'k');
    ht = plot(8, 0);                     % dummy handle
    %set(gca, 'XTick', []); set(gca, 'YTick', []);
    xlabel(''); ylabel('');
    title('Control the snake with arrow keys. Quit with any other key.');
    axis([-1, 1, -1, 1, 0, maxf]); shg;
    t = 1;                               % convex factor for nodes
    tic;
    go = plot(.7*scribble('ready?'), 'r', LW, lw);
    shg
    pause(1)
    delete(go);
    go = plot(.4*scribble('go!'), 'r', LW, lw);
    shg
    pause(.5)
    delete(go);

    while ( d ~= 0 )                     % until game over or quit
        
        t = t + .51*alfa;
        if ( t > 1 )
            t = 0; dr = res*d;
            s = [ s(2:end), s(end)+dr ];
            if ( length(s) < len + pts )
                s = [ s(2), s ]; 
            end
        end
        
        y = (1-t)*s(1:end-1)+t*s(2:end);
        if ( nodes == 1 )
            c = chebfun(y.');
            c = c(chebpts(5*length(y)));
        elseif ( nodes == 2 )
            fhd = min(ceil(0.4*sqrt(length(y))), 4);
            c = bary(linspace(-1, 1, 5*length(y)), y, linspace(-1, 1, ...
                length(y)), weights(length(y)-1, fhd));
        elseif ( nodes == 0 )
            yy = linspace(-1, 1, length(y)).';
            c = polyfit(yy, y, length(y) - 1, dom);
        end
        
        for k = 1:numel(hs1)
            delete(hs1(k));
        end
        
        hs1 = plot3(real(c), imag(c), f(real(c), imag(c)), 'b-', LW, lw);
        delete(hs2);
        delete(hs1s);
        delete(hs2s);
        hs2 = [plot3(real(y(1:end-1)), imag(y(1:end-1)), f(real(y(1:end-1)), ...
            imag(y(1:end-1))), 'bo', 'MarkerFaceColor', 'b', LW, lw), ...
            plot3(real(y(end)), imag(y(end)), f(real(y(end)), imag(y(end))), ...
            'bo', LW, lw)];
        hs1s = plot(c, 'k-', LW, lw/2);
        hs2s = plot(y, 'k.', LW, lw/2);
        drawnow;
        
        tic
        % check if the snake hits itself or the boundary
        if ( max(abs([real(y(end)), imag(y(end))])) > 1 || ...
                min(abs(y(end)-y(1:end-1))) < res/2 )
            ht = plot(.8*scribble('game over'), 'r', LW, lw);
            chebtune(dd, .5);
            shg; pause(1);
            fails = fails + 1;
            if ( fails > failmax )
                d = 0; 
            end
            break
        end
        
        if ( abs(y(end)-fd) < res/2 ) % snake eats food ?
            pts = pts + 1;
            chebtune(grb, .5);
            if ( ~rem(pts, 10) )
                lvl = lvl + 1;
                alfa = alfa * 1.1;
                chebtune(10*chebpoly(pts));
            end
            if ( ~rem(pts, 30) && (pts - 1) )
                fails = fails - 1;
                up = plot(.8*scribble('1 up!'), 'r', LW, lw);
                chebtune(lv, 1);
                shg; pause(1);
                delete(up);
            end
            title(['Points : ' num2str(pts) '       Level : ' num2str(lvl) ...
                '       Lives: ' num2str(failmax-fails)], 'color', 'k');
            fd = food();
            while ( any( abs(fd-y) < res/2) )
                fd = food();
            end
            set(hf, 'XData', real(fd), 'YData', imag(fd), 'ZData', f(real(fd), imag(fd)));
            set(hfs, 'XData', real(fd), 'YData', imag(fd));
            if ~rem(pts+1, 30) && (pts-1)
                set(hf, 'MarkerFaceColor', 'r', 'Color', [0, 0.6, 0]);
            else
                set(hf, 'MarkerFaceColor', [0, 0.6, 0], 'Color', 'r');
            end
        end
        
    end
    
    for k = 1:numel(ht)
        delete(ht(k));
    end
    
end

plot(.8*scribble('goodbye'), 'r', LW, lw); chebtune(kld, 1);
shg; pause(1); close(gcf);
warning(W);

    function w = weights(n, fhd) % weights for Floater-Hormann interpolation
        w = zeros(1, n+1);
        for l = 0:n
            ji = max(l-fhd, 0);
            jf = min(l, n-fhd);
            sumcoeff = zeros(jf-ji+1, 1);
            for i=ji:jf
                sumcoeff(i-ji+1) = nchoosek(fhd, l-i);
            end
            w(l+1) = (-1)^(l-fhd)*sum(sumcoeff);
        end
    end

end
