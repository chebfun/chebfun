function [k,v]=findpeaks2(x,m,w)
%FINDPEAKS finds peaks with optional quadratic interpolation [K,V]=(X,M,W)
%
%  Inputs:  X        is the input signal (does not work with UInt datatype)
%           M        is mode:
%                       'q' performs quadratic interpolation
%                       'v' finds valleys instead of peaks
%           W        is the width tolerance; a peak will be eliminated if there is
%                    a higher peak within +-W samples
%
% Outputs:  K        are the peak locations in X (fractional if M='q')
%           V        are the peak amplitudes: if M='q' the amplitudes will be interpolated
%                    whereas if M~='q' then V=X(K).

% Outputs are column vectors regardless of whether X is row or column.
% If there is a plateau rather than a sharp peak, the routine will place the
% peak in the centre of the plateau. When the W input argument is specified,
% the routine will eliminate the lower of any pair of peaks whose separation
% is <=W; if the peaks have exactly the same height, the second one will be eliminated.
% All peak locations satisfy 1<K<length(X).
%
% If no output arguments are specified, the results will be plotted.
%
%       Copyright (C) Mike Brookes 2005
%      Version: $Id: findpeaks.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    m=' ';
end
nx=length(x);
if any(m=='v')
    x=-x(:);        % invert x if searching for valleys
else
        x=x(:);        % force to be a column vector
end
dx=x(2:end)-x(1:end-1);
r=find(dx>0);
f=find(dx<0);

if length(r)>0 & length(f)>0    % we must have at least one rise and one fall
    dr=r;
    dr(2:end)=r(2:end)-r(1:end-1);
    rc=repmat(1,nx,1);
    rc(r+1)=1-dr;
    rc(1)=0;
    rs=cumsum(rc); % = time since the last rise
    df=f;
    df(2:end)=f(2:end)-f(1:end-1);
    fc=repmat(1,nx,1);
    fc(f+1)=1-df;
    fc(1)=0;
    fs=cumsum(fc); % = time since the last fall
    rp=repmat(-1,nx,1);
    rp([1; r+1])=[dr-1; nx-r(end)-1];
    rq=cumsum(rp);  % = time to the next rise
    fp=repmat(-1,nx,1);
    fp([1; f+1])=[df-1; nx-f(end)-1];
    fq=cumsum(fp); % = time to the next fall
    k=find((rs<fs) & (fq<rq) & (floor((fq-rs)/2)==0));   % the final term centres peaks within a plateau
    v=x(k);
    if any(m=='q')         % do quadratic interpolation
        b=0.5*(x(k+1)-x(k-1));
        a=x(k)-b-x(k-1);
        j=(a>0);            % j=0 on a plateau
        v(j)=x(k(j))+0.25*b(j).^2./a(j);
        k(j)=k(j)+0.5*b(j)./a(j);
        k(~j)=k(~j)+(fq(k(~j))-rs(k(~j)))/2;    % add 0.5 to k if plateau has an even width
    end
    % now purge nearby peaks
    
    if nargin>2
        j=find(k(2:end)-k(1:end-1)<=w);
        while any(j)
            j=j+(v(j)>=v(j+1));
            k(j)=[];
            v(j)=[];
            j=find(k(2:end)-k(1:end-1)<=w);
        end
    end
else
    k=[];
    v=[];
end
if any(m=='v')
    v=-v;    % invert peaks if searching for valleys
end
if ~nargout
    if any(m=='v')
        x=-x;    % re-invert x if searching for valleys
        ch='v';
    else
        ch='^';
    end
    plot(1:nx,x,'-',k,v,ch);
end