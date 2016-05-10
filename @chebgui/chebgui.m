classdef chebgui
% INTRODUCTION TO CHEBGUI
%
% Chebgui is a graphical user interface to Chebfun's capabilities for solving
% ODEs and PDEs (ordinary and partial differential equations) and eigenvalue
% problems. More precisely, it deals with the following classes of problems.  In
% all cases both single equations and systems of equations can be treated, as
% well as integral and integro-differential equations.
%
% ODE BVP (boundary-value problem): an ODE or system of ODEs on an interval
% [a,b] with boundary conditions at both boundaries.
%
% ODE IVP (initial-value problem): an ODE or system of ODEs on an interval [a,b]
% with boundary conditions at just one boundary. By default, IVPs are
% automatically reformulated to coupled first order systems, which are then
% solved using MATLAB's built in solver ODE113.
%
% ODE EIGENVALUE PROBLEM: a differential or integral operator or system of
% operators on an interval [a,b] with homogeneous boundary conditions, where we
% want to compute one or more eigenvalues and eigenfunctions. Generalized
% eigenvalue problems L*u = lambda*M*u are also treated.
%
% PDE BVP: a time-dependent problem of the form u_t = N(u,x,t), where N is a
% linear or nonlinear differential operator.
%
% For ODEs, Chebgui assumes that the independent variable, which varies over the
% interval [a,b], is either r, t or x, and that the dependent variable(s) have
% name(s) different from r, t and x.
%
% For eigenvalue problems, Chebgui assumes that the eigenvalue is called l, lam
% or lambda.
%
% For PDEs, Chebgui assumes that the space variable on [a,b] is x and that the
% time variable, which ranges over a span t1:dt:t2 is t.
%
% Here is a three-sentence sketch of how the solutions are computed.  The ODE
% and eigenvalue problems are solved by Chebfun's automated Chebyshev spectral
% methods underlying the Chebfun commands <backslash> and SOLVEBVP.  The
% discretizations involved will be described in a forthcoming paper by Driscoll
% and Hale, and the Newton and damped Newton iterations used to handle
% nonlinearity is described in [1]. The PDE problems are solved by Chebfun's
% PDE15S method, due to Hale, which is based on spectral discretization in x
% coupled with Matlab's ODE15S solution in t.
%
% To use Chebgui, the simplest method is to type chebgui at the Matlab prompt.
% The GUI will appear with a demo already loaded and ready to run; you can get
% its solution by pressing the green SOLVE button.  To try other preloaded
% examples, open the DEMOS menu at the top.  To input your own example on the
% screen, change the windows in ways which we hope are obvious. Inputs are
% vectorized, so x*u and x.*u are equivalent, for example.  Derivatives are
% indicated by single or double primes, so a second derivative is u'' or u".
%
% The GUI allows various types of syntax for describing the differential
% equation and boundary conditions of problems. The differential equations can
% either be in anonymous function form, e.g.
%
%   @(u) diff(u,2)+x.*sin(u)
%
% or a "natural syntax form", e.g.
%
%   u''+x.*sin(u)
%
% The first format gives extra flexibility, e.g. for specifying an integral
% operator with the help of CUMSUM.
%
% Boundary conditions can be in either of these two forms, or alternatively one
% can specify homogeneous Dirichlet or Neumann conditions simply by typing
% 'dirichlet' or 'neumann' in the relevant fields.  Eigenvalue problems can be
% specified by equations like
%
%   -u" + x^2*u = lambda*u
%
% The GUI supports systems of coupled equations, where the problem can most
% easily be set up with a format like
%
%   u' + sin(v) = u+v
%   cos(u) + v' = 0
%
% or
%
%   u" = lambda*v
%   v" = lambda*(u+u')
%
% It is possible to right-click on an input box, which brings up a larger input
% box to make it easier to input more complicated problem. Observe that on Apple
% computers, CTRL+Click is equivalent to right-clicking.
%
%
% Finally, the most valuable Chebgui capability of all is the "Export to m-file"
% button.  With this feature, you can turn an ODE or PDE solution from the GUI
% into an M-file in standard Chebfun syntax.  This is a great starting point for
% more serious explorations.
%
% CHEBGUI is also the constructor for chebgui objects. For example
%    chebg = chebgui('type','bvp','domain','[-1,1]', ...
%                    'de','u" = sin(u)','lbc','u = 1','rbc','u = 0')
%    show(chebg)
%
% See also chebop/solvebvp, chebop/eigs, chebfun/pde15s, chebfun/ode45,
% chebfun/ode113, chebfun/ode15s, chebfun/bvp4c, chebfun/bvp5c.
%
% References:
%   [1] A. Birkisson and T. A. Driscoll, “Automatic Fréchet Differentiation for
%   the Numerical Solution of Boundary-Value Problems,” ACM Transactions on
%   Mathematical Software, vol. 38, no. 4, Article 26, Aug. 2012.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        type = '';          % Type of chebgui (bvp, pde, or eig)
        domain = '';        % Spacial domain (may contain breakpoints)
        DE = '';            % Differential equation, or rhs in u_t = ... for PDEs
        LBC = '';           % Left BCs
        RBC = '';           % Right BCs
        BC = '';            % Interior / nonstandard BCs
        timedomain = '';    % Time domain (should include breakpoints)
        sigma = '';         % Third input to an EIGS call
        init = '';          % Initial guess/condition for nonlin BVPs/PDEs
        tol = '5e-13';      % Solution tolerance

        % This initalises the OPTIONS struct for CHEBGUI. It containts a list of
        % miscellaneous options for when solving problems with CHEBGUI, namely:
        %   damping:        Whether the Newton iteration is to be damped or not.
        %   grid:           Show grids on plots in CHEBGUI.
        %   discretization: Discretization for ODEs ('values' or 
        %                   'coeffs').
        %   pdeholdplot:    Hold plot during solving PDEs.
        %   fixYaxisLower:  Fix lower y axis while solving PDEs.
        %   fixYaxisUpper:  Fix upper y axis while solving PDEs.
        %   fixN:           Fix discretization for when solving PDEs.
        %   numeigs:        Number of eigenvalues sought when EIG problems are
        %                   solved.
        options = struct('damping', '1', ...
            'plotting', '0.5', ...
            'grid', 1, ...
            'discretization', 'values', ...
            'ivpSolver', 'ode113', ...
            'pdeSolver', 'pde15s', ...
            'pdeholdplot', 0, ...
            'fixYaxisLower', '',  ...
            'fixYaxisUpper', '', ...
            'fixN', '', ...
            'numeigs', '');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function c = chebgui(varargin)
        %CHEBGUI    The CHEBGUI constructor
        %
        % The command CHEBGUI in the command window opens the CHEBGUI
        % figure window, with a random BVP example shown.
        %
        % One can also pass arguments to the CHEBGUI constructor to create
        % CHEBGUI objects. Here, the call will be of the form
        %
        %   CG = CHEBGUI('field', value)
        %
        % To see the possible fields of the CHEBGUI objects, call
        %   DISP(CHEBGUI)
            
            % No input --> load random example to the GUI window
            if ( isempty(varargin) )
                c = chebgui.demo();
                show(c);    % Open the GUI
                return      % No further action needed
            end
            
            % Look at the first input:
            v1 = varargin{1};
            
            % Single input --> load existing chebgui object/file
            if ( nargin == 1 )
                if ( isa(v1, 'chebgui') )
                    % Calling chebgui constructor with a chebgui loads the GUI.
                    c = v1;
                elseif ( ischar(v1) )
                    % Calling with a single string loads a .guifile`
                    if ( ~exist(v1,'file') )
                        % Throw a method if this .guifie is not found
                        error('CHEBFUN:CHEBGUI:chebgui:missingFile', ...
                              'Unable to find file: %s', v1);
                    end
                    % Load the existing demp
                    c = chebgui.demo2chebgui(v1);
                end
                if ( nargout == 0 )
                    show(c); % Open the GUI
                end
                return
            end
            
            % Multiple inputs --> loop through them (using CHEBGUI/SET)
            k = 1;
            while ( k < nargin )
                c = set(c, varargin{k:k+1});
                k = k + 2;
            end

        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    methods ( Access = public, Static = false )

        function allVarNames = getVarNames(guifile)
            %GETVARNAMES    Return a cell-array with all variables names
            %
            % ALLVARNAMES = GETVARNAMES(GUIFILE) returns a cell-array, whose
            % entries correspond to the names of all variables that appear in
            % a problem specified by the CHEBGUI object GUIFILE.
            
            % Look at the input for the differential equation:
            deInput = guifile.DE;
            
            % Wrap the input strings in a cell (if it's not a cell already)
            if ( isa(deInput, 'char') )
                deInput = cellstr(deInput);
            end
            
            % Obtain useful strings describing the differential equation part:
            [dummy, dummy, dummy, dummy, dummy, dummy, allVarNames] = ...
                setupFields(guifile, deInput, 'DE');
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Convert initial conditions to a useful format
        bcOut = bcReform(dom, bcInput, isIorF)
        
        % Information shown during BVP solving.
        [dummy, displayTimer] = displayBVPinfo(handles, mode, varargin)
        
        % Information shown after IVP solving.
        [dummy, displayTimer] = displayIVPinfo(handles, u, isIVP, varargin)
        
        % Return a random BVP CHEBGUI demo.
        cg = demo()

        % Load a demo stored in a .guifile to a CHEBGUI object
        cg = demo2chebgui(demoPath)
        
        % Construct initial guesses/conditions prescribed by CHEBGUI objects
        init = constructInit(initInput, allVarNames, indVarNameSpace, xt);
    end
         
end
