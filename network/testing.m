
m = 4;
x = linspace(-10,10,m);
y = x;
z = x;

[X,Y,Z] = ndgrid(x,y,z);

constants.epsilon = 1/4/pi;
constants.alpha = 4*pi;  % solid angle for 0D nodes is 4 pi since no fluid is blocked by the geometry
constants.mu = 1; % viscosity of Newtonian fluid

Stokeslets.n_nodes = numel(X);
Stokeslets.nodes = [X(:) Y(:) Z(:)];

% links = NaN(3*m^2*(m-1),2); % number of links from Wrobel et al 2016
Stokeslets.links = [];
sz = size(X);
for s = 1:Stokeslets.n_nodes
    [i,j,k] = ind2sub(sz,s);
    
    stencil = [i-1, j, k;  i+1,j,k;  i,j-1,k;  i,j+1,k;  i,j,k-1;  i,j,k+1;];
    for st = 1:size(stencil,1)
        if any(stencil(st,:) < 1) || any(stencil(st,:) > sz) % stencil point doesn't lie within range of actual grid
            continue
        end
        Stokeslets.links(end+1,:) = [sub2ind(sz,stencil(st,1),stencil(st,2),stencil(st,3)) s];
        
    end
end

Stokeslets.links = unique(sort(Stokeslets.links,2),'rows');
Stokeslets.n_links = size(Stokeslets.links,1);
Stokeslets.E = repmat(1,Stokeslets.n_links,1);
% Stokeslets.l_0 = sqrt(sum( (Stokeslets.nodes(Stokeslets.links(:,1),:)  - Stokeslets.nodes(Stokeslets.links(:,2),:) ).^2 , 2 ) );
Stokeslets.l_0 =  sqrt(sum(  (( Stokeslets.nodes(Stokeslets.links(:,1),:)  + Stokeslets.nodes(Stokeslets.links(:,2),:) ) / 2).^2 , 2));
% resting lengths are the initial distances from the center of each link to the origin
Stokeslets.l = Stokeslets.l_0;
Stokeslets.eta = repmat(10,Stokeslets.n_links,1);

for i = 1:Stokeslets.n_nodes
    Stokeslets.link_members{i,1} = [];
    for j = 1:Stokeslets.n_links
        [ism,ind] = ismember(i,Stokeslets.links(j,:));
        if ism
            Stokeslets.link_members{i,1}(1,end+1) = j; % link index that node i is a member of
            if ind == 1
                Stokeslets.link_members{i,1}(2,end) = 1; % if node is first member of link, force on this node = f_s
            else
                Stokeslets.link_members{i,1}(2,end) = -1; % if node is 2nd member of link, force on this node = -f_s
            end
        end
    end
end

%%

d = Stokeslets.nodes(Stokeslets.links(:,2),:)  - Stokeslets.nodes(Stokeslets.links(:,1),:) ;
r = sqrt(sum( d.^2 , 2 ) );

f_s = Stokeslets.l_0.^2 .* Stokeslets.E .* (r./Stokeslets.l - 1) .* d ./ r;


g = NaN(Stokeslets.n_nodes,3); % net Maxwell element viscoelastic force on each network node due to all attached links
u = NaN(Stokeslets.n_nodes,3);

parfor i = 1:Stokeslets.n_nodes
    g(i,:) = sum(   Stokeslets.link_members{i}(2,:)' .*  f_s( Stokeslets.link_members{i}(1,:) , :) , 1);
    
end

parfor i = 1:Stokeslets.n_nodes
    temp = 0;
    for j = 1:Stokeslets.n_nodes
        temp = temp + g(j,:) *  calcS(Stokeslets.nodes(j,:)',Stokeslets.nodes(i,:)',epsilon^2);
    end
    u(i,:) = 1/mu * 1/alpha * 1/2 * temp;
end


dldt = Stokeslets.E .* Stokeslets.l_0 ./ Stokeslets.eta .* ( r ./ Stokeslets.l - 1);


figure(124);
plot3(Stokeslets.nodes(:,1),Stokeslets.nodes(:,2),Stokeslets.nodes(:,3),'ko','MarkerFaceColor','k');
text(Stokeslets.nodes(:,1)+0.2,Stokeslets.nodes(:,2),Stokeslets.nodes(:,3),cellfun(@num2str,num2cell(1:Stokeslets.n_nodes),'UniformOutput',false),'FontSize',12);
hold on
plot3([Stokeslets.nodes(Stokeslets.links(:,1),1) Stokeslets.nodes(Stokeslets.links(:,2),1)]',...
    [Stokeslets.nodes(Stokeslets.links(:,1),2) Stokeslets.nodes(Stokeslets.links(:,2),2)]',...
    [Stokeslets.nodes(Stokeslets.links(:,1),3) Stokeslets.nodes(Stokeslets.links(:,2),3)]',...
    'r--','LineWidth',2);

quiver3(Stokeslets.nodes(:,1),Stokeslets.nodes(:,2),Stokeslets.nodes(:,3),u(:,1),u(:,2),u(:,3),'b-','LineWidth',2.5);
quiver3(Stokeslets.nodes(:,1),Stokeslets.nodes(:,2),Stokeslets.nodes(:,3),g(:,1),g(:,2),g(:,3),'b--','LineWidth',2.5);

hold off

xlabel('x'); ylabel('y'); zlabel('z');
axis equal
return
%%

fun = @(t,y) derivatives(t,y,Stokeslets,constants);

y0 = [ Stokeslets.l; reshape(Stokeslets.nodes',Stokeslets.n_nodes*3,1) ];
tspan = [0 0.1];

[t,y] = ode45(fun, tspan, y0,odeset('Refine',8));
return
%%
clear y
t = [0:0.1:10];
y(1,:) = y0;

for i = 2:length(t)
    
y(i,:) = y(i-1,:)' + (t(i) - t(i-1)) * fun(t(i),y(i-1,:)');
end


%%
clear x
l = y(:,1:Stokeslets.n_links);
for i = 1:size(y,1)
% x(:,:,i) = reshape( y(i,Stokeslets.n_links + 1 : Stokeslets.n_links + Stokeslets.n_nodes*3) , 3,[])';
x = reshape( y(i,Stokeslets.n_links + 1 : Stokeslets.n_links + Stokeslets.n_nodes*3) , 3,[])';

figure(123);
plot3(x(:,1),x(:,2),x(:,3),'ko','MarkerFaceColor','k');
% text(Stokeslets.nodes(:,1)+0.2,Stokeslets.nodes(:,2),Stokeslets.nodes(:,3),cellfun(@num2str,num2cell(1:Stokeslets.n_nodes),'UniformOutput',false),'FontSize',12);
hold on
plot3([x(Stokeslets.links(:,1),1) x(Stokeslets.links(:,2),1)]',...
    [x(Stokeslets.links(:,1),2) x(Stokeslets.links(:,2),2)]',...
    [x(Stokeslets.links(:,1),3) x(Stokeslets.links(:,2),3)]',...
    'r--','LineWidth',2);

hold off

xlabel('x'); ylabel('y'); zlabel('z');
axis equal
xlim([-25 25]);  ylim([-25 25]);  zlim([-25 25]); 
title(num2str(t(i)));
% pause(0.1)
drawnow

end