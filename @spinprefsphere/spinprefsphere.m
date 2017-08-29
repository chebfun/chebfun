classdef spinprefsphere < spinpreference
%SPINPRESPHERE   Class for managing preferences when solving a PDE on the
%sphere with SPINSPHERE.
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
%   grid                      * If 'on', plots longitude and latitude circles.                      
%     ['off']                    
%      'on'
%
%   iterplot                  * Plots the solution every ITERPLOT iterations of
%     [1]                       the time-stepping loop if 'plot' is 'movie'.
%
%
%   Nplot                     * Number of grid points in each direction for 
%     [256]                     the movie. If Nplot>N, the data are interpolated 
%                               to a finer grid.
%
%   plot                      * Plot options: 'movie' to plot a movie of the
%     ['movie']                 the solution, 'off' otherwise. 
%      'off'
%
%   scheme                    * Time-stepping scheme. HELP/IMEX for the list
%     []                        of available schemes. Default is empty, i.e., 
%                               automatically chosen by the code depending on
%                               the PDE.
%
%   view                      * Viewpoint specification when 'plot' is 'movie'.
%     [-37.5 30]   
%
% Construction:
%
%   PREF = SPINPREFSPHERE() creates a SPINPREF2 object with the default values.
%
%   PREF = SPINPREFSPHERE(PROP1, VALUE1, PROP2, VALUE2, ...) creates a 
%   SPINPREFSPHERE object with the properties PROP1 and PROP2 set to VALUE1 and
%   VALUE2.
%
% See also SPINSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        Clim                  % Limits of the colorbar (1x2*NVARS DOUBLE)
        colormap              % Color look-up table (STRING)
        grid                  % For longitude and latitude circles (STRING)
        view = [-37.5 30];    % Viewpoint of the plot (1x2 DOUBLE)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinprefsphere(varargin) 
            if ( nargin == 0 )
                pref.colormap = 'parula';
                pref.grid = 'off';
                pref.dataplot = 'real';
                pref.dealias = 'off';
                pref.iterplot = 1;
                pref.Nplot = 128;
                pref.plot = 'movie';
            else
                pref = spinprefsphere();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end