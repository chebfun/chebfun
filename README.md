About
=====

Chebfun is a collection of algorithms and an open-source software system in
object-oriented MATLAB which extends familiar powerful methods of numerical
computation involving numbers to continuous or piecewise-continuous functions.
It also implements continuous analogues of linear algebra notions like the QR
decomposition and the SVD, and solves ordinary differential equations. The
mathematical basis of the system combines tools of Chebyshev expansions, fast
Fourier transform, barycentric interpolation, recursive zerofinding, and
automatic differentiation.


Installation and requirements
=============================

Chebfun Version 5.0 is compatible with MATLAB 7.6 (2008a) and above.

Installation instructions:

1. Unzip the contents of the zip file to a directory. We suggest the name
   "chebfun" for this directory. Make sure to maintain the existing
   subdirectory structure of the zip package. (Please note: If you install
   into the "toolbox" subdirectory of the MATLAB program hierarchy, you will
   need to click the button "Update toolbox path cache" from the
   File/Preferences... dialog in MATLAB.)

2. In MATLAB, add the chebfun directory to your path. This can be done by
   selecting *File/Set Path...* from the main or command window menus, or with
   the command `pathtool`. We recommend that you select the *Save* button on
   this dialog so that Chebfun is on the path automatically in future MATLAB
   sessions. (Alternatively, you can put an `addpath` command in your
   `startup.m` file, as described in the MATLAB documentation.)

3. A restart of MATLAB may be needed if you want to access the user guides via
   the Help browser, or a call to `clear classes` if you had a previous
   version of Chebfun installed.

Point your web browser to the guide/html directory of the Chebfun installation
in order to read the user guides.

Please see the file LICENSE.txt for licensing information.

Chebfun can be found on the web at http://www.chebfun.org/


Getting started
===============

[todo]


Contributing
============

[todo]


License
=======

See `LICENSE.txt` for Chebfun's licensing informtion.
