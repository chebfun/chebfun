function varargout = surf(F, varargin)
%SURF Plot of the surface represented by a CHEBFUN2V.
%
% SURF(F) is the surface plot of F, where F is a CHEBFUN2V with three
% components.
%
% SURF(F,'-') also shows the seams of the parameterisation on the surface.
%
% SURF(F,...) allows for the same plotting options as Matlab's SURF
% command.
%
% See also CHEBFUN2/SURF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 100; 

% Empty check:
if ( isempty(F) )
    surf([]);
    return
end

ish = ishold;

nF = F.nComponents;

if ( nF == 2 )
    error('CHEBFUN:CHEBFUN2V:surf:notEnoughComponents', ...
        'CHEBFUN2V does not represent a surface as it has only two components');
end

% Making varargin non-empty
if ( isempty(varargin) )
    varargin = {};
end

% plot seams in coordinate representation.
if ( ~isempty(varargin) )
    if ( length(varargin{1})<5 )
        % Only option with <=3 letters is a colour, marker, line
        ll = regexp(varargin{1},'[-:.]+','match');
        cc = regexp(varargin{1},'[bgrcmykw]','match');  % color
        
        if ( ~isempty(ll) )
            if ( strcmpi(ll{1},'.') )
                % We have . first. Don't plot a line.
                ll = {};
            elseif ( strcmpi(ll{1},'.-') )
                % We have . first. Don't plot a line.
                ll{1} = '-';
            end
        end
        plotline = ~isempty(ll);  % plot row and col pivot lines?
        
        % default to black seam.
        if ( isempty(cc) ) 
            cc{1} = 'k';
        end
        
        % call to CHEBFUN2/SURF
        F1 = F.components{1};
        F2 = F.components{2};
        F3 = F.components{3};
        h1 = surf( F1, F2, F3, varargin{2:end} ); hold on

        if ( plotline )
            LW = 'linewidth'; lw = 2;
            dom = F1.domain;
            
            x = chebpts( numpts, dom(1:2) );
            lft = dom(3) * ones(length(x), 1);
            h2 = plot3(F1(x, lft), F2(x ,lft), F3(x,lft), ...
                   'linestyle', ll{1}, 'Color', cc{1}, LW, lw ); hold on
            
            rght = dom(4) * ones(length(x), 1);
            h3 = plot3(F1(x, rght), F2(x, rght), F3(x, rght),...
                    'linestyle', ll{1}, 'Color', cc{1}, LW, lw );
            
            y = chebpts(numpts, dom(3:4) );
            dwn = dom(1) * ones( length(x), 1 );
            h4 = plot3(F1(dwn, y), F2(dwn, y), F3(dwn, y),...
                   'linestyle', ll{1}, 'Color', cc{1}, LW, lw );
            
            up = dom(2) * ones(length(x), 1);
            h5 = plot3(F1(up, y), F2(up, y), F3(up, y),...
                   'linestyle', ll{1}, 'Color', cc{1}, LW, lw );
            
            h = [h1 h2 h3 h4 h5];
        else
            h = h1;
        end
    else
        % straight call to CHEBFUN2/SURF
        h = surf(F.components{1}, F.components{2}, F.components{3}, varargin{:});
    end
else
    % straight call to CHEBFUN2/SURF
    h = surf(F.components{1}, F.components{2}, F.components{3}, varargin{:});
end

if ( ~ish )
    hold off
end

if ( nargout > 1 )
    varargout = {h};
end

end
