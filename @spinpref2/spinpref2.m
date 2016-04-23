classdef spinpref2 < spinpreference
%SPINPRE2   Class for managing preferences when solving a 2D PDE with SPIN2.
%
% Available preferences ([] = defaults):
%
%   Clim                      * Limits of the colorbar when 'plot' is 'movie'.
%     []                        Default is empty, i.e., automatically chosen by 
%                               the code. 
%
%   dataToPlot                * What data to plot when the solution is complex-
%     ['real']                  valued.
%      'imag'
%      'abs'
%
%   dealias                   * If it is 'on', use the 2/3-rule to zero high 
%     ['off']                   wavenumbers.
%      'on'
% 
%   dt                        * Time-step for time discretization. Default is 
%     []                        empty, i.e., adaptively chosen by the code to 
%                               achieve errTol. 
%
%   dtmax                     * Maximum time-step when using an apative grid in
%     [5]                       time.
%
%   dtmin                     * Minimum time-step when using an apative grid in
%     [1e-10]                   time.
%
%   errTol                    * Desired accuracy on the solution.
%     [1e-4]
%
%   iterPlot                  * Plot the solution every ITERPLOT iterations of
%     [1]                       the time-stepping loop if 'plot' is 'movie'.
%
%   M                         * Number of points for complex means to evaluate
%     [32]                      the phi-functions.
%
%   N                         * Number points in each direction for spatial 
%     [64]                      discretization. To switch to adaptive grid, set
%                               N=[].
%
%   Nmin                      * Minimum number of points in each direction when 
%     [64]                      using an adaptive grid in space.
%
%   Nmax                      * Maximum number of points in each direction when   
%     [512]                     using an adaptive grid in space.
%                                    
%   Nplot                     * Number of grid points in each direction for 
%     [256]                     plotting. If Nplot>N, the data are interpolated 
%                               to a finer grid.
%
%   plot                      * Plot options: 'movie' to plot a movie of the
%     ['movie']                 the solution, 'off' otherwise. 
%      'off'
%
%   scheme                    * Time-stepping scheme. HELP/SPINPSCHEME for the
%     ['etdrk4']                list of available schemes.
%
%   view                      * Viewpoint specification when 'plot' is 'movie'.
%     [0 90]   
%
% Construction:
%
%   PREF = SPINPREF2() creates a SPINPREF2 object with the default values.
%
%   PREF = SPINPREF2(PDECHAR) creates a SPINPREF2 object corresponding to the 
%   preferences used for the SPIN2(PDECHAR) demo. Strings available include
%   'GL2' for Ginzburg-Landau equation and 'GS2' for Gray-Scott equations. 
%   Other PDEs are available, see HELP/SPIN2.
%
%   PREF = SPINPREF2(PROP1, VALUE1, PROP2, VALUE2, ...) creates a SPINPREF2 
%   object with the properties PROP1 and PROP2 set to VALUE1 and VALUE2.
%
% See also SPIN2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        Clim                  % Limits of the colorbar (1x2*NVARS DOUBLE)
        view = [0 90];        % Viewpoint of the plot (1x2 DOUBLE)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref2(varargin) 
            if ( nargin == 0 )
                pref.dataToPlot = 'real';
                pref.dealias = 'off';
                pref.dtmin = 1e-10;
                pref.dtmax = 5;
                pref.errTol = 1e-4;
                pref.iterPlot = 1;
                pref.M = 32;
                pref.N = 64;
                pref.Nmin = 64;
                pref.Nmax = 512;
                pref.Nplot = 256;
                pref.plot = 'movie';
                pref.scheme = 'etdrk4';
            elseif ( nargin == 1 )
                pdechar = varargin{1};
                pref.dataToPlot = 'real';
                pref.dealias = 'off';
                pref.dtmin = [];
                pref.dtmax = [];
                pref.errTol = 1e-4;
                pref.M = 32;
                pref.Nmin = [];
                pref.Nmax = [];
                pref.Nplot = 256;
                pref.plot = 'movie';
                pref.scheme = 'etdrk4';
                if ( strcmpi(pdechar, 'GL2') == 1 )
                    pref.dt = 2e-1;
                    pref.iterPlot = 1;
                    pref.N = 64;
                elseif ( strcmpi(pdechar, 'GS2') == 1 )
                    pref.dt = 4;
                    pref.iterPlot = 8;
                    pref.N = 64;
                elseif ( strcmpi(pdechar, 'Schnak2') == 1 )
                    pref.dt = 5e-1;
                    pref.iterPlot = 10;
                    pref.N = 64;
                elseif ( strcmpi(pdechar, 'SH2') == 1 )
                    pref.dt = 1;
                    pref.iterPlot = 1;
                    pref.N = 64; 
                else
                    error('SPINPREF2:CONSTRUCTOR', 'Unrecognized PDE.')
                end
            else
                pref = spinpref2();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end