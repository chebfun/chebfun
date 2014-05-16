classdef chebgui

% TODO:  Reread this documentation and make sure it makes sense in v5.

% INTRODUCTION TO CHEBGUI
%
% Chebgui is a graphical user interface to Chebfun's capabilities for
% solving ODEs and PDEs (ordinary and partial differential equations) and
% eigenvalue problems. More precisely, it deals with the following classes
% of problems.  In all cases both single equations and systems of equations
% can be treated, as well as integral and integro-differential equations.
%
% ODE BVP (boundary-value problem): an ODE or system of ODEs on an interval
% [a,b] with boundary conditions at both boundaries.
%
% ODE IVP (initial-value problem): an ODE or system of ODEs on an interval
% [a,b] with boundary conditions at just one boundary. (However, for
% complicated IVPs like the Lorenz equations, other methods such as
% chebfun/ode45 will be much more effective.)
%
% ODE EIGENVALUE PROBLEM: a differential or integral operator or system of
% operators on an interval [a,b] with homogeneous boundary conditions,
% where we want to compute one or more eigenvalues and eigenfunctions.
% Generalized eigenvalue problems L*u = lambda*M*u are also treated.
%
% PDE BVP: a time-dependent problem of the form u_t = N(u,x,t), where N is
% a linear or nonlinear differential operator.
%
% For ODEs, Chebgui assumes that the independent variable, which varies
% over the interval [a,b], is either x or t, and that the dependent
% variable(s) have name(s) different from x and t.
%
% For eigenvalue problems, Chebgui assumes that the eigenvalue is called
% lambda or lam or l.
%
% For PDEs, Chebgui assumes that the space variable on [a,b] is x and that
% the time variable, which ranges over a span t1:dt:t2 is t.
%
% Here is a three-sentence sketch of how the solutions are computed.  The
% ODE and eigenvalue problems are solved by Chebfun's automated Chebyshev
% spectral methods underlying the Chebfun commands <backslash> and
% SOLVEBVP.  The discretizations involved will be described in a
% forthcoming paper by Driscoll and Hale, and the Newton and damped Newton
% iterations used to handle nonlinearity will be described in a forthcoming
% paper by Birkisson and Driscoll.  The PDE problems are solved by
% Chebfun's PDE15S method, due to Hale, which is based on spectral
% discretization in x coupled with Matlab's ODE15S solution in t.
%
% To use Chebgui, the simplest method is to type chebgui at the Matlab
% prompt.  The GUI will appear with a demo already loaded and ready to run;
% you can get its solution by pressing the green SOLVE button.  To try
% other preloaded examples, open the DEMOS menu at the top.  To input your
% own example on the screen, change the windows in ways which we hope are
% obvious. Inputs are vectorized, so x*u and x.*u are equivalent, for
% example.  Derivatives are indicated by single or double primes, so a
% second derivative is u'' or u".
%
% The GUI allows various types of syntax for describing the differential
% equation and boundary conditions of problems. The differential equations
% can either be in anonymous function form, e.g.
%
%   @(u) diff(u,2)+x.*sin(u)
%
% or a "natural syntax form", e.g.
%
%   u''+x.*sin(u)
%
% The first format gives extra flexibility, e.g. for specifying an
% integral operator with the help of CUMSUM.
%
% Boundary conditions can be in either of these two forms, or alternatively
% one can specify homogeneous Dirichlet or Neumann conditions simply by
% typing 'dirichlet' or 'neumann' in the relevant fields.  Eigenvalue
% problems can be specified by equations like
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
% Finally, the most valuable Chebgui capability of all is the "Export to
% m-file" button.  With this feature, you can turn an ODE or PDE solution
% from the GUI into an M-file in standard Chebfun syntax.  This is a great
% starting point for more serious explorations.
%
% CHEBGUI is also the constructor for chebgui objects. For example
%    chebg = chebgui('type','bvp','domain','[-1,1]', ...
%                    'de','u" = sin(u)','lbc','u = 1','rbc','u = 0')
%    show(chebg)
%
% See also chebop/solvebvp, chebop/eigs, chebfun/pde15s, chebfun/ode45,
% chebfun/ode113, chebfun/ode15s, chebfun/bvp4c, chebfun/bvp5c.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    properties
        type = '';          % Type of chebgui (bvp, pde, or eig)
        domain = '';        % Spacial domain (may contain breakpoints)
        DE = '';            % Differential equation, or rhs in u_t = ... for PDEs
        LBC = '';           % Left BCs
        RBC = '';           % Right BCs
        BC = '';            % Interior / nonstandard BCs
        timedomain = '';    % Time domain (should include breakpoints)
        sigma = '';         % Third input to an EIGS call
        init = '';          % Initial guess/condition for nonlin BVPs/PDEs
        tol = 1e-10;        % Solution tolerance

	% TODO:  Documentation.
        options = struct('damping', '1', ...
	    'plotting', '0.5', ...
	    'grid', 1, ...
            'discretization', @colloc2, ...
            'pdeholdplot', 0, ...
	    'fixYaxisLower', '',  ...
            'fixYaxisUpper', '', ...
	    'fixN', '', ...
	    'numeigs', '');
    end
    
    methods
        function c = chebgui(varargin)

	% TODO:  Document inputs to this more carefully.
            
            % No input --> load random example to the GUI window
            if ( isempty(varargin) )
                c = loadexample(chebgui('type', 'bvp'), -1);
                show(c); % Open the GUI
                return
            end
            
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
                        error('CHEBFUN:chebgui:missingfile', ...
                              'Unable to find file: %s', v1);
                    end
                    % Load the existing demp
                    c = loaddemos(chebgui('type', 'bvp'), v1);
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
    
    methods (Static = true)
        
        [dummy, displayTimer] = displayBVPinfo(handles, mode, varargin);
        
        
        function initialiseFigures(handles)
            %INITIALISEFIGURES      Reset figures in the CHEBGUI window.

            % Reset the top figures
            cla(handles.fig_sol, 'reset');
            title('Solutions')
	    box on

            % Reset the bottom figure
            cla(handles.fig_norm, 'reset');
            title('Updates')
	    box on
        end
        
    end
    
    methods( Access = private )
       
        % Export a BVP to an .m file
        exportBVP2mfile(guifile, pathname, filename)
        
        % Export an EIG problem to an .m file
        exportEIG2mfile(guifile, pathname, filename, handles)
        
        % Export a PDE to an .m file
        exportPDE2mfile(guifile, pathname, filename)
        
        % Convert a string to an anonymous function
        varargout = strConv2anon(guifile, str, type)
        
        % Load the CHEBGUI demos
        cg = loaddemos(guifile, guifilepath)
        
        % Populate the fields of the GUI
        [field, allVarString, indVarName, pdeVarNames, pdeflag, eigVarNames, ...
            allVarNames]  = setupFields(guifile, input, type, allVarString)
        
        % Solve a GUI BVP
        varargout = solveGUIbvp(guifile,handles)
        
        % Solve a GUI EIG problem
        varargout = solveGUIeig(guifile,handles)
    end
    
    methods( Access = private, Static = true )
        
        % The CHEBGUI lexer
        [out, varNames, pdeVarNames, eigVarNames, indVarNames] = ...
            strConvLexer(str, type)
        
        % The LL(1) parser
        parseOut = strConvParser(lexIn)
        
        % Convert a syntax tree to prefix format
        prefixOut = strConvTree2prefix(syntaxTree)
        
        % Get rid of unecessary parenthesis in infix format strings
        strOut = strConvParSimp(strIn)
        
        % Split expressions separated by commas to individual syntax trees
        treeOut = splitTree_commas_equalSigns(treeIn)
        
        % Split a BVP syntax tree
        treeOut = splitTree_bvp(treeIn)
        
        % Split an EIG syntax tree
        [newTree, lambdaTree, lambdaSign] = splitTree_eig(treeIn)
        
        % Split a PDE syntax tree
        [newTree, pdeSign] = splitTree_pde(treeIn)
        
        % Convert an expression on PREFIX form to INFIX form
        [infixOut, notaVAR] = strConvPref2inf(prefIn)
        
    end
    
end
