
var = 'U';
ind = 2;


x = [Solutions.phase; 2*pi]';
 y = [Solutions.(var)(:,ind); Solutions.(var)(1,ind)]' ;
%  y(3) = NaN;
 
 x2 = linspace(0,2*pi,1000);
 
 
 curve2 = csape(x,y,'periodic');
 
 tempx = x(1:end-1);  tempy = y(1:end-1);
%  tempx(3) = [];  tempy(3) = [];
 
 tfit = trig_interp_fit(tempx, tempy);
 
%  curve1 = trig_interp_eval(interpolant(3).(var)(ind),x2);
 curve1 = trig_interp_eval(tfit,x2);
 
 curve1 = [curve1 , curve1];
 xplt = [x2 , x2 + 2*pi];
 
 curve2 = fnval(curve2,x2);
 curve2 = [curve2, curve2];
 
 
 figure(384)
 plot(x,y,'o','markerfacecolor','b');
  hold on
 plot(x + 2*pi , y , 'bo');

 plot(xplt,curve1,'b-');
 plot(xplt,curve2,'r-');
 hold off
 grid on