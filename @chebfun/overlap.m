function [f, g] = overlap(f, g)
%OVERLAP  Overlap the domain of two chebfun objects.
%
%   [FOUT, GOUT] = OVERLAP(F ,G) returns two chebfuns such that FOUT.domain ==
%   GOUT.domain and F(x) = FOUT(x), G(x) = GOUT(x) for all x in domain of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab some data:
fDom = f.domain; 
gDom = g.domain;
fImps = f.impulses;
gImps = g.impulses;
fRows = size(fImps, 1); 
gRows = size(gImps, 1);
maxRows = max(fRows, gRows);

if ( length(fDom) == length(gDom) && all(fDom == gDom) )
    % Domains already match. Just make impulses long enough.
    f.impulses = [ fImps ; zeros(maxRows - fRows, length(fDom))];
    g.impulses = [ gImps ; zeros(maxRows - gRows, length(fDom))];    

else
    %% Check domain is valid:
    hs = max(f.hscale, g.hscale);
    if ( norm([fDom(1) - gDom(1), fDom(end) - gDom(end)], inf) > 1e-15*hs )
        error('CHEBFUN:overlap:domains', ...
            'Inconsistent domains; domain(f) ~= domain(g).')
    end
    % Take the union of the domains:
    newDom = union(fDom, gDom);
    
    %% F
    % Make the new funs for f:
    fOut = cell(1, numel(newDom) - 1 );
    c = 1; % counter
    for k = 1:numel(f.funs) % Loop through existing funs.
        % Find breaks to introduce in this interval:
        dk = f.funs{k}.domain;
        newdk = newDom(newDom >= dk(1) & newDom <= dk(2));
        % Call restrict on the funs to output an array:
        fOut(1,c:c+numel(newdk)-2) = restrict(f.funs{k}, newdk);
        % Increment the counter:
        c = c + numel(newdk) - 1;
    end
    
    %% G
    % Make the new funs for g:
    gOut = cell(1, numel(newDom) - 1 );
    c = 1; % counter
    for k = 1:numel(g.funs) % Loop through existing funs.
        % Find breaks to introduce in this interval:
        dk = g.funs{k}.domain;
        newdk = newDom(newDom >= dk(1) & newDom <= dk(2));
        % Call restrict on the funs to output an array:
        gOut(1,c:c+numel(newdk)-2) = restrict(g.funs{k}, newdk);
        % Increment the counter:
        c = c + numel(newdk) - 1;
    end

    %% Impulses!
    fOutImps = zeros(maxRows, length(newDom));
    [ignored, fIndex, fOutInd] = intersect(fDom, newDom);
    
    fOutImps(2:fRows, fOutInd) = fImps(2:fRows, fIndex);
    idx = abs(fOutImps(2:end,:)) > 100*eps;
    idx = sum(idx, 1) ~= 0; % Indices with deltas.
    % Compute the average for delta indices:
    if ( any(idx) )
        fOutImps(1,idx) = (feval(f, newDom(idx), 'left') + ...
                           feval(f,newDom(idx),'right')) / 2;
    end
    % Otherwise compute the normal function values:
    if ( ~all(idx) )
        fOutImps(1,~idx) = feval(f, newDom(~idx));
    end
    
    gOutImps = zeros(maxRows, length(newDom)); 
    [ignored, gIndex, gOutInd] = intersect(gDom, newDom);
    gOutImps(2:gRows,gOutInd) = gImps(2:gRows,gIndex);
    idx = abs(gOutImps(2:end,:)) > 100*eps;
    idx = sum(idx,1) ~= 0; % Indices with deltas.
    % Compute the average for delta indices:
    if ( any(idx) )
        gOutImps(1,idx) = (feval(g,newDom(idx), 'left') + ...
                           feval(g,newDom(idx),'right')) / 2;
    end
    % Otherwise compute the normal function values:
    if ( ~all(idx) )
        gOutImps(1,~idx) = feval(g, newDom(~idx));
    end
    
    %% Append new funs/domains/impulses to f and g for output:
    f.funs = fOut; 
    f.domain = newDom; 
    f.impulses = fOutImps;
    g.funs = gOut; 
    g.domain = newDom; 
    g.impulses = gOutImps;

end
