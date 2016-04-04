function [err, time] = spincomp(pdechar, tspan, u0, pref)
%SPINCOMP  Compare time-stepping schemes.
%   SPINCOMP(PDECHAR, TSPAN, U0, PREF) solves the PDE specified by the STRING
%   PDECHAR on TSPAN x U0.DOMAIN, with initial condition U0, using the various
%   time-steps DT stored in PREF.DT and the time-stepping schemes listed in
%   PREF.SCHEME. PREF is a SPINPREF OBJECT. See HELP/SPIN, HELP/SPIN2 and 
%   HELP/SPIN3 for a list of available strings.
%
% Example 1: Compare the ETDRK schemes for 1D Kuramoto-Sivashinsky equation
%
%    dom = [0 32*pi]; tspan = [0 30];
%    u0 = chebfun('cos(x/16).*(1 + sin(x/16))',dom,'trig');
%    pref = spinpref();
%    pref.dt = [1, 5e-1, 2e-1, 1e-1, 5e-2, 2e-2, 1e-2];
%    pref.scheme = {'etdrk4', 'exprk5s8', 'krogstad', 'friedli', ...
%       'hochbruck-ostermann', 'minchev', 'strehmel-weiner'};
%    pref.N = 256;
%    spincomp('KS', tspan, u0, pref)
%
% Example 2: Compare the PREDICTOR-CORRECTOR schemes for 2D Gray-Scott equations
%
%    G = 1.25; dom = G*[0 1 0 1]; tspan = [0 30];
%    u01 = chebfun2(@(x,y) 1-exp(-150*((x-G/2).^2 + (y-G/2).^2)), dom, 'trig');
%    u02 = chebfun2(@(x,y) exp(-150*((x-G/2).^2 + 2*(y-G/2).^2)), dom, 'trig');
%    u0 = chebmatrix(u01);
%    u0(2,1) = u02;
%    pref = spinpref2;
%    pref.dt = [1, 5e-1, 2e-1, 1e-1, 5e-2, 2e-2];
%    pref.scheme = {'pec423', 'pecec433', 'pec524', 'pecec534', 'pec625', ...
%       'pecec635', 'pec726', 'pecec736'};
%    pref.N = 128;
%    spincomp('GS2', tspan, u0, pref)
%
% See also SPIN, SPIN2, SPIN3, SPINPREF, SPINPREF2, SPINPREF3, SPINSCHEME.
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the schemes:
schemes = pref.scheme;
nSchemes = length(schemes);

% Get the time-steps:
timesteps = pref.dt;
nTimesteps = length(timesteps);

% Convert the initial conditiona to a CHEBMATRIX:
if ( isa(u0, 'chebfun') == 1 || isa(u0, 'chebfun2') == 1 || ...
    isa(u0, 'chebfun3') == 1)
    u0 = chebmatrix(u0);
elseif ( isa(u0, 'chebfun2v') == 1 || isa(u0, 'chebfun3v') == 1 )
    temp = chebmatrix(u0(1));
    for k = 2:size(u0, 1)
        temp(k,1) = u0(k);
    end
    u0 = temp;
end
nVars = length(u0);

% Create a grid to compare exact and computed solutions:
N = 32;
dom = u0{1}.domain;
xx = trigpts(N, dom);
if ( isa(u0{1}, 'chebfun') == 1 )
    dim = 1;
elseif ( isa(u0{1}, 'chebfun2') == 1 )
    dim = 2;
    [xx, yy] = meshgrid(xx);
elseif ( isa(u0{1}, 'chebfun3') == 1 )
    dim = 3;
    [xx, yy, zz] = meshgrid(xx);
end
    
% First, estimate the exact solution using a very small time step (half the 
% smallest time-step) and with PECEC736 (7th-order multistep scheme):
if ( isa(pref, 'spinpref') == 1 )
    prefuexact = spinpref();
elseif ( isa(pref, 'spinpref2') == 1 )
    prefuexact = spinpref2();
elseif ( isa(pref, 'spinpref3') == 1 )
    prefuexact = spinpref3();
end
prefuexact.scheme = 'pecec736';
prefuexact.dt = min(timesteps)/2;
prefuexact.plot = 'off';
prefuexact.M = pref.M;
prefuexact.N = pref.N;
uexact = spinoperator.solvepde(pdechar, tspan, u0, prefuexact);
if ( isa(uexact, 'chebmatrix') == 0 )
    uexact = chebmatrix(uexact);
end

% Scale (i.e., maximum amplitude) of the exact solution:
scale = 0;
for i = 1:nVars
    if ( dim == 1 )
        scale = max(scale, max(abs(uexact{i}(xx))));
    elseif( dim == 2 )
        scale = max(scale, max(max(abs(uexact{i}(xx,yy)))));
    elseif ( dim == 3 )
        scale = max(scale, max(max(max(abs(uexact{i}(xx,yy,zz))))));
    end
end

% Second, compute solutions for different schemes and time-steps:
prefu = spinpref();
prefu.plot = 'off';
prefu.M = pref.M;
prefu.N = pref.N;
err = zeros(nTimesteps, nSchemes);
time = zeros(nTimesteps, nSchemes);
for k = 1:nTimesteps
    for l = 1:nSchemes
        prefu.scheme = schemes{l};
        prefu.dt = timesteps(k);
        tic, u = spinoperator.solvepde(pdechar, tspan, u0, prefu);  
        time(k,l) = toc;
        if ( isa(u, 'chebmatrix') == 0 )
            u = chebmatrix(u);
        end
        isNan = isNanTest(u, dim);
        if ( isNan == 1 )
            err(k,l) = NaN;
        else
            for i = 1:nVars
                if ( dim == 1 )
                    temp = abs(u{i}(xx) - uexact{i}(xx));
                elseif( dim == 2 )
                    temp = abs(u{i}(xx,yy) - uexact{i}(xx,yy));
                elseif ( dim == 3 )
                    temp = abs(u{i}(xx,yy,zz) - uexact{i}(xx,yy,zz));
                end
                temp = max(temp(:));
                err(k,l) = max(err(k,l), temp);
            end
        end
    end
end
err = err/scale; % Relative error
TF = tspan(end);
timesteps = timesteps/TF; % Relative time-steps

% Plot Accuarcy vs time-step:
figure, subplot(1, 2, 1)
labels = cell(nSchemes, 1);
for l = 1:nSchemes
    labels{l} = schemes{l};
    if ( l > 1 )
        hold on
    end
    if ( l < 8 )
        loglog(timesteps, err(:, l), '.-', 'linewidth', 2, 'markersize', 30)
    elseif ( l < 15 )
        loglog(timesteps, err(:, l), 'x-', 'linewidth', 2, 'markersize', 10)
    else
        loglog(timesteps, err(:, l), 'o-', 'linewidth', 2, 'markersize', 10)
    end
end
left = 10^floor(log10(min(timesteps)));
right = 10^ceil(log10(max(timesteps)));
down = max(min(1e-10, 10^floor(log10(min(err(:))))), 1e-12);
up = min(max(1e0, 10^ceil(log10(max(err(:))))), 1e2);
set(gca, 'fontsize', 16), axis([left right down up])
xlabel('Relative time-step'), ylabel(sprintf('Relative error at t = %.3f', TF))
legend(labels, 'Location', 'NorthWest')

% Plot Accuarcy vs computer time:
subplot(1, 2, 2)
for l = 1:nSchemes
    if ( l > 1 )
        hold on
    end
    if ( l < 8 )
        loglog(time(:,l), err(:, l), '.-', 'linewidth', 2, 'markersize', 30)
   elseif ( l < 15 )
        loglog(time(:,l), err(:, l), 'x-', 'linewidth', 2, 'markersize', 10)
    else
        loglog(time(:,l), err(:, l), 'o-', 'linewidth', 2, 'markersize', 10)
    end
end
left = 10^(floor(log10(min(min(time)))));
right = 10^(ceil(log10(max(max(time)))));
set(gca, 'fontsize', 16), axis([left right down up])
xlabel('Computer time (s)'), ylabel(sprintf('Relative error at t = %.3f', TF))
legend(labels, 'Location', 'NorthEast')

end

function out = isNanTest(u, dim)

B = u.blocks;
if ( dim == 1 )
    out = any(cell2mat(cellfun(@(C) isnan(C), B, 'UniformOutput', 0))); 
elseif ( dim == 2 )
    out = any(cell2mat(cellfun(@(C) isnan(C.cols), B, 'UniformOutput', 0))); 
    out = out || ...
        any(cell2mat(cellfun(@(C) isnan(C.rows), B, 'UniformOutput', 0))); 
elseif ( dim == 3 )
    out = any(cell2mat(cellfun(@(C) isnan(C.cols), B, 'UniformOutput', 0))); 
    out = out || ...
        any(cell2mat(cellfun(@(C) isnan(C.rows), B, 'UniformOutput', 0))); 
    out = out || ...
        any(cell2mat(cellfun(@(C) isnan(C.tubes), B, 'UniformOutput', 0))); 
end

end
