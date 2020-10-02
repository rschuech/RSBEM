function hh = confellipse2(xy,conf)
%CONFELLIPSE2 Draws a confidence ellipse.
% CONFELLIPSE2(XY,CONF) draws a confidence ellipse on the current axes
% which is calculated from the n-by-2 matrix XY and encloses the
% fraction CONF (e.g., 0.95 for a 95% confidence ellipse).
% H = CONFELLIPSE2(...) returns a handle to the line.

% written by Douglas M. Schwarz
% schwarz@kodak.com
% last modified: 12 June 1998

n = size(xy,1);
mxy = mean(xy);

numPts = 181; % The number of points in the ellipse.
th = linspace(0,2*pi,numPts)';


p = 2; % Dimensionality of the data, 2-D in this case.

k = finv(conf,p,n-p)*p*(n-1)/(n-p);
% Comment out line above and uncomment line below to use ftest toolbox.
% k = fdistinv(p,n-p,1-conf)*p*(n-1)/(n-p);

[pc,score,lat] = princomp(xy);
% Comment out line above and uncomment 3 lines below to use ftest toolbox.
% xyp = (xy - repmat(mxy,n,1))/sqrt(n - 1);
% [u,lat,pc] = svd(xyp,0);
% lat = diag(lat).^2;

ab = diag(sqrt(k*lat));
exy = [cos(th),sin(th)]*ab*pc' + repmat(mxy,numPts,1);

% Add ellipse to current plot
h = line(exy(:,1),exy(:,2),'Clipping','off');
if nargout > 0
    hh = h;
end