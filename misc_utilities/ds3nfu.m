function  [Xf, Yf] = ds3nfu(HAX, X, Y, Z) ;
% Takes X,Y, and Z coordinates on the axis HAX, and gives you the
% X and Y position in the figure window, from 0 to 1
% same idea as ds2nfu.m from file exchange, except it works with 3d input
% -Jeff Marks

%%% Get the relative positions in X and Y axis coordinates.
[Xr,Yr,Zr] = ds3nfu_relative(HAX,X,Y,Z)

% find the current units, set to normalized, then set back to the current
% units
Units = get(HAX,'Units');
set(HAX,'Units','normalized')
P = get(HAX,'Position');
set(HAX,'Units',Units)

% get the minimum and maximum limits of the axis in figure units
Xmin = P(1);
Xmax = P(1)+P(3);
Ymin = P(2);
Ymax = P(2)+P(4);

% transform the axis corners into relative figure space
V = axis(HAX);
if length(V) == 6
    Xcorn = V([1 1 1 1 2 2 2 2]);
    Ycorn = V([3 3 4 4 3 3 4 4]);
    Zcorn = V([5 6 5 6 5 6 5 6]);
else 
    Xcorn = V([1 1 1 1 2 2 2 2]);
    Ycorn = V([3 3 4 4 3 3 4 4]);
    Zcorn = zeros(1,8);
end
% transform the corners of the axes into relative figure units
[Xc,Yc,Zc] = ds3nfu_relative(HAX,Xcorn,Ycorn,Zcorn);
% find the X and Y limits
Xcmin = min(Xc);
Xcmax = max(Xc);
Ycmin = min(Yc);
Ycmax = max(Yc);

Xscale = (Xmax-Xmin)/(Xcmax-Xcmin);
Yscale = (Ymax-Ymin)/(Ycmax-Ycmin);

% transform the relative figure units into absolute figure units
Xf = Xmin + Xscale*(Xr-Xcmin);
Yf = Ymin + Yscale*(Yr-Ycmin);

function [Xrel,Yrel,Zrel] = ds3nfu_relative(HAX,X,Y,Z)
V = axis(HAX);
[AZ,EL] = view(HAX);

% rotation is valid only if centered on axes and scaled equally
% relative to the axis.  It is okay to go ahead and
% transform, because these are all relative units.
if length(V)==6 & abs(EL)~=90
    Center = [mean(V(1:2)) mean(V(3:4)) mean(V(5:6))];
    Scale = [V(2)-V(1)   V(4)-V(3)  V(6)-V(5)];
    X = (X - Center(1))./Scale(1);
    Y = (Y - Center(2))./Scale(2);
    Z = (Z - Center(3))./Scale(3);
else
 
    Center = [mean(V(1:2)) mean(V(3:4))];
    Scale = [V(2)-V(1)   V(4)-V(3) ];
    % if EL is 90 or -90, AZ behaves differently and must be scaled itself
    % tan does not preserve information if abs AZ is >90
    if abs(AZ)>90
        overtan = 1;
    else
       overtan = 0;
    end
    AZ = AZ*pi/180;
    YTAN = tan(AZ);
    if overtan
        XTAN = -1;
        YTAN = -YTAN; 
    else
        XTAN = 1;
    end
    YTAN2 = YTAN * Scale(2)/Scale(1);
    AZ = 180* atan2(YTAN2,XTAN)/pi;

    X = (X - Center(1))./Scale(1);
    Y = (Y - Center(2))./Scale(2);
    Z = zeros(size(X));
end


if EL==-90 & ~any(abs(get(HAX,'CameraUpVector'))==1)
    % Matlab flips the axes only on these conditions
    X = -X;
    Y = -Y;
end

[TH,R] = cart2pol(X,Y);
TH2 = TH - AZ*pi/180;
[X4,Y2] = pol2cart(TH2,R);

[TH3,R3] = cart2pol(Y2,Z);
TH4 = TH3 + EL*pi/180;
[Y4,Z4] = pol2cart(TH4,R3);

Xrel = X4;
%%% The Y and Z axis need to be swiched, because AZ=0 and EL=0 have Z as
%%% the up-down axis and X as the left-right axis
Yrel = Z4;
Zrel = Y4;
