classdef spinpref < spinpreference
%SPINPREF   Class for managing preferences when solving a 1D PDE with SPIN.
%
% Available preferences ([] = defaults):
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
%     [20]                      the time-stepping loop if 'plot' is 'movie'.
%
%   M                         * Number of points for complex means to evaluate
%     [64]                      the phi-functions.
%       
%   Nplot                     * Number of grid points for the movie. If Nplot>N,
%     [1024]                    the data are interpolated on a finer grid.
%    
%   plot                      * Plot options: 'movie' to plot a movie of the 
%     ['movie']                 solution, 'waterfall' to use the WATERFALL
%      'waterfall'              command, 'off' otherwise.
%      'off'
%
%   scheme                    * Time-stepping scheme. HELP/EXPINTEG for the
%     ['etdrk4']                list of available schemes.
%
%   Ylim                      * Limits of the y-axis when 'plot' is 'movie'.
%     []                        Default is empty, i.e., automatically chosen by 
%                               the code. 
%              
% Construction:
%
%   PREF = SPINPREF() creates a SPINPREF object with the default values.
%
%   PREF = SPINPREF(PROP1, VALUE1, PROP2, VALUE2, ...) creates a SPINPREF object
%   with the properties PROP1 and PROP2 set to VALUE1 and VALUE2.
%
% See also SPIN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        M            % Number of points for complex means (1x1 INT)
        Ylim         % Limit of the y-axis of the plot (1x2*NVARS DOUBLE)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref(varargin) 
            if ( nargin == 0 )
                pref.dataplot = 'real';
                pref.dealias = 'off';
                pref.iterplot = 20;
                pref.M = 64;
                pref.Nplot = 1024;
                pref.plot = 'movie';
                pref.scheme = 'etdrk4';
            else
                pref = spinpref();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end