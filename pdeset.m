function varargout = pdeset(varargin)
%PDESET    Set options for PDE15S
%   PDESET('NAME1', VALUE1, 'NAME2', VALUE2,...) creates options for the
%   CHEBFUN/PDE15S() routine. It acts as a gateway to ODESET() for the usual ODE
%   options for use in advancing through time, in addition to some new options.
%
%   OPTIONS = PDESET(OLDOPTS, 'NAME1', VALUE1,...) alters an existing options
%   structure OLDOPTS.
%
%   PDESET PROPERTIES (In addition to ODESET properties)
%
%       Eps - Tolerance to use in solving the PDE [ positive scalar {1e-6} ].
%
%       N - Turn off spacial adaptivity. [{NaN} | positive integer  ]
%           Use a fixed spacial grid of size N. If N is NaN, then the automatic
%           procedure is used.
%
%       Plot - Plot the solution at the end of every time chunk. [ {on} | off ]
%              Turning this off can improve speed considerably.
%
%       HoldPlot - Hold the plots after each chunk. [ on | {off} ]
%
%       YLim - Fix the limits of the Y axis if plotting. [ 2x1 vector | {NaN} ]
%           If Ylim is NaN then the imits are determined automatically.
%   
%       PlotStyle - Change the plotting options. [ string | ''-'' ].
%
%       PDEflag - Specify which entries correspond to time derivatives. 
%           [ vector of logicals {true} ].
%
%       AdjustBCs - Adjust boundary conditions to be satisfied by initial
%                   conditions.
%           [ logical {true} ].
%
% See also ODESET.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%  PDE specific options
names = ['Eps      ' 
         'N        '
         'Plot     '
         'HoldPlot '
         'YLim     '
         'PlotStyle'
         'PDEflag  '
         'AdjustBCs']; 
     
m = size(names, 1);
shortNames = cell(m, 1);
for k = 1:m
    shortNames{k} = strtrim(names(k, :));
end

% Initialise:
opts = {};
pdeOpts = {};

if ( nargin == 0 )
    if ( nargout == 0 )
        odeset;
        fprintf('             Eps: [ positive scalar {1e-6} ]\n')
        fprintf('               N: [ {NaN} | positive integer  ]\n')        
        fprintf('            Plot: [ {on} | off ]\n')
        fprintf('        HoldPlot: [ on | {off} ]\n')
        fprintf('            YLim: [ 2x1 vector | {NaN} ]\n')
        fprintf('       PlotStyle: [ string | ''-'']\n')
        fprintf('         PDEflag: [ vector of logicals {true} ]\n')
        fprintf('       AdjustBCs: [ logical {true} ]\n')
    else
        % Get the ODE opts:
        opts = odeset();
        % Add empty PDE opts:
        for j = 1:m
            opts.(shortNames{j}) = [];
        end
        varargout{1} = opts;
    end      

    return
end

% Is an ODESET() / PDESET() structure being passed?
if ( isstruct(varargin{1}) )
    opts = varargin{1};
    varargin(1) = [];
end

% Remember the old PDEOPT values:
for k = 1:m
    namek = shortNames{k};
    if ( isfield(opts, namek) )
        pdeOpts = [ pdeOpts {namek, opts.(namek)}];
    end
end

% Parse the remaining input and update pdeopts entries
k = 1;
while ( k < length(varargin) )
    if ( ~any(strcmpi(fieldnames(odeset), varargin{k})) )
        if ( strcmpi(varargin{k}, 'Plot') || strcmpi(varargin{k}, 'HoldPlot') )
            varargin{k + 1} = onoff(varargin{k + 1});
        end
        pdeOpts = [pdeOpts varargin(k : k + 1)];
        varargin(k : k + 1) = [];
    else
        k = k + 2;
    end
end

% Get the ODE opts:
opts = odeset(opts, varargin{:});

% Add empty PDE opts:
for j = 1:m
    opts.(shortNames{j}) = [];
end

% Attach the PDE opts:
for k = 1:2:length(pdeOpts)
    for j = 1:m
        if ( strcmpi(pdeOpts{k}, shortNames{j}) )
            opts.(shortNames{j}) = pdeOpts{k + 1};
            break
        end
        if ( j == m )
            error('CHEBFUN:pdeset:unknownOption', ...
                ['Unrecognized property name ', pdeOpts{k}, '.'])
        end
    end
end

% Output the opts:
varargout{1} = opts;

end

function out = onoff(in)
% Convert logical values to 'on' or 'off'.
    if ( ~ischar(in) )
        if ( logical(in) )
            out = 'on';
        else
            out = 'off';
        end
    else
        out = in;
    end
end
