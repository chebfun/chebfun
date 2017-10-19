classdef spinpref3 < spinpreference
%SPINPREF3   Class for managing preferences when solving a 3D PDE with SPIN3.
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
%     [128]                     the movie. If Nplot>N, the data are interpolated 
%                               on a finer grid.
%
%   plot                      * Plot options: 'movie' to plot a movie of the
%     ['movie']                 solution, 'off' otherwise.
%      'off'
%
%   scheme                    * Time-stepping scheme. HELP/EXPINTEG for the
%     ['etdrk4']                list of available schemes.
%
%   slices                    * Slices of the volumetric slice plot when 'plot'
%     []                        is 'movie'. Default is empty, i.e., 
%                               automatically chosen by the code.
%                
% Construction:
%
%   PREF = SPINPREF3() creates a SPINPREF3 object with the default values.
%
%   PREF = SPINPREF3(PROP1, VALUE1, PROP2, VALUE2, ...) creates a SPINPREF3
%   object with the properties PROP1 and PROP2 set to VALUE1 and VALUE2.
%
% See also SPIN3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        Clim         % Limits of the colorbar (1x2*NVARS DOUBLE)
        colormap     % Color look-up table (STRING)
        M            % Number of points for complex means (1x1 INT)
        slices       % Slices of the volumetric slice plot (1x3 CELL, slices{1}
                     % is a DOUBLE of positions x corresponding to the slices 
                     % x=cst, same for slices{2} and slices{3} with y and z)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref3(varargin) 
            if ( nargin == 0 )
                pref.colormap = 'parula';
                pref.dataplot = 'real';
                pref.dealias = 'off';
                pref.iterplot = 1;
                pref.M = 32;
                pref.Nplot = 64;
                pref.plot = 'movie';
                pref.scheme = 'etdrk4';
            else
                pref = spinpref3();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end