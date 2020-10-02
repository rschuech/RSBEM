


% load temporal S/N Results file here

body = [6 0.5]; %body of interest to examine tail guesses for

bodies = [ [Results.AR1]'  [Results.AR2]'];

inds = ismember(bodies, body,'rows');


tails = [ [Results.amp]' [Results.lambda]'  [Results.nlambda]' ];
tails = tails(inds,:);
obj = [Results.taxis];  obj = [obj.temporal];
for o = 1:numel(obj)
    if isempty(obj(o).SN)
        obj(o).SN = NaN;
    end
end
obj = [obj.SN]';
obj = obj(inds);


F = scatteredInterpolant(tails(:,1),tails(:,2),obj,'linear','none');
npts = 200;
 [X,Y] = meshgrid(linspace(min(tails(:,1)),max(tails(:,1)),npts),linspace(min(tails(:,2)),max(tails(:,2)),npts));
 Z = F(X,Y);
 
 figure(640)
 pcolor(X,Y,Z);  shading interp
 hold on
 dots = plot(tails(:,1),tails(:,2),'ko','markerfacecolor','k','markersize',8);
 hold off
 grid on
 colorbar