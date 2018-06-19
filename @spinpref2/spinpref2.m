classdef spinpref2 < spinpreference
%SPINPRE2   Class for managing preferences when solving a 2D PDE with SPIN2.
%
% Available preferences ([] = defaults):
%
%   Clim                      * Limits of the colorbar when 'plot' is 'movie'.
%     []                        Default is empty, i.e., automatically chosen by 
%                               the code. 
%
%   colormap                  * Color look-up table. See HELP/COLORMAP.
%     ['parula'] 
%
%   dataplot                  * Plotting options when the solution is complex.
%     ['real']                 
%      'imag'
%      'abs'
%
%   dealias                   * If 'on', uses the 2/3-rule to zero out high
%     ['off']                   wavenumbers.
%      'on'
%
%   iterplot                  * Plots the solution every ITERPLOT iterations of
%     [1]                       the time-stepping loop if 'plot' is 'movie'.
%
%   M                         * Number of points for complex means to evaluate
%     [32]                      the phi-functions.
%
%   Nplot                     * Number of grid points in each direction for 
%     [256]                     the movie. If Nplot>N, the data are interpolated 
%                               on a finer grid.
%
%   plot                      * Plot options: 'movie' to plot a movie of the
%     ['movie']                 the solution, 'off' otherwise. 
%      'off'
%
%   scheme                    * Time-stepping scheme. HELP/EXPINTEG for the
%     ['etdrk4']                list of available schemes.
%
%   view                      * Viewpoint specification when 'plot' is 'movie'.
%     [0 90]   
%
% Construction:
%
%   PREF = SPINPREF2() creates a SPINPREF2 object with the default values.
%
%   PREF = SPINPREF2(PROP1, VALUE1, PROP2, VALUE2, ...) creates a SPINPREF2 
%   object with the properties PROP1 and PROP2 set to VALUE1 and VALUE2.
%
% See also SPIN2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        Clim                  % Limits of the colorbar (1x2*NVARS DOUBLE)
        colormap              % Color look-up table (STRING)
        M                     % Number of points for complex means (1x1 INT)
        view = [0 90];        % Viewpoint of the plot (1x2 DOUBLE)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref2(varargin) 
            if ( nargin == 0 )
                pref.colormap = 'parula';
                pref.dataplot = 'real';
                pref.dealias = 'off';
                pref.iterplot = 1;
                pref.M = 32;
                pref.Nplot = 128;
                pref.plot = 'movie';
                pref.scheme = 'etdrk4';
            else
                pref = spinpref2();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end