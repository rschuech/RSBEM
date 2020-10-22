
m = 4;
x = linspace(-5,5,m);
y = x;
z = x;

[X,Y,Z] = ndgrid(x,y,z);

constants.epsilon = 1/4/pi;
constants.alpha = 4*pi;  % solid angle for 0D nodes is 4 pi since no fluid is blocked by the geometry
constants.mu = 1; % viscosity of Newtonian fluid

Network.n_nodes = numel(X);
Network.nodes = [X(:) Y(:) Z(:)];

% links = NaN(3*m^2*(m-1),2); % number of links from Wrobel et al 2016
Network.links = [];
sz = size(X);
for s = 1:Network.n_nodes
    [i,j,k] = ind2sub(sz,s);
    
    stencil = [i-1, j, k;  i+1,j,k;  i,j-1,k;  i,j+1,k;  i,j,k-1;  i,j,k+1;];
    for st = 1:size(stencil,1)
        if any(stencil(st,:) < 1) || any(stencil(st,:) > sz) % stencil point doesn't lie within range of actual grid
            continue
        end
        Network.links(end+1,:) = [sub2ind(sz,stencil(st,1),stencil(st,2),stencil(st,3)) s];
        
    end
end

Network.links = unique(sort(Network.links,2),'rows');
Network.n_links = size(Network.links,1);
Network.E = repmat(1,Network.n_links,1);
% Network.l_0 = sqrt(sum( (Network.nodes(Network.links(:,1),:)  - Network.nodes(Network.links(:,2),:) ).^2 , 2 ) );
Network.l_0 =  sqrt(sum(  (( Network.nodes(Network.links(:,1),:)  + Network.nodes(Network.links(:,2),:) ) / 2).^2 , 2));
% resting lengths are the initial distances from the center of each link to the origin
Network.l = Network.l_0;
Network.eta = repmat(10,Network.n_links,1);

for i = 1:Network.n_nodes
    Network.link_members{i,1} = [];
    for j = 1:Network.n_links
        [ism,ind] = ismember(i,Network.links(j,:));
        if ism
            Network.link_members{i,1}(1,end+1) = j; % link index that node i is a member of
            if ind == 1
                Network.link_members{i,1}(2,end) = 1; % if node is first member of link, force on this node = f_s
            else
                Network.link_members{i,1}(2,end) = -1; % if node is 2nd member of link, force on this node = -f_s
            end
        end
    end
end

%%

d = Network.nodes(Network.links(:,2),:)  - Network.nodes(Network.links(:,1),:) ;
r = sqrt(sum( d.^2 , 2 ) );

f_s = Network.l_0.^2 .* Network.E .* (r./Network.l - 1) .* d ./ r;


g = NaN(Network.n_nodes,3); % net Maxwell element viscoelastic force on each network node due to all attached links
u = NaN(Network.n_nodes,3);

parfor i = 1:Network.n_nodes
    g(i,:) = sum(   Network.link_members{i}(2,:)' .*  f_s( Network.link_members{i}(1,:) , :) , 1);
    
end

parfor i = 1:Network.n_nodes
    temp = 0;
    for j = 1:Network.n_nodes
        temp = temp + g(j,:) *  calcS(Network.nodes(j,:)',Network.nodes(i,:)',epsilon^2);
    end
    u(i,:) = 1/mu * 1/alpha * 1/2 * temp;
end


dldt = Network.E .* Network.l_0 ./ Network.eta .* ( r ./ Network.l - 1);


figure(124);
plot3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),'ko','MarkerFaceColor','k');
text(Network.nodes(:,1)+0.2,Network.nodes(:,2),Network.nodes(:,3),cellfun(@num2str,num2cell(1:Network.n_nodes),'UniformOutput',false),'FontSize',12);
hold on
plot3([Network.nodes(Network.links(:,1),1) Network.nodes(Network.links(:,2),1)]',...
    [Network.nodes(Network.links(:,1),2) Network.nodes(Network.links(:,2),2)]',...
    [Network.nodes(Network.links(:,1),3) Network.nodes(Network.links(:,2),3)]',...
    'r--','LineWidth',2);

% quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.u(:,1),Network.u(:,2),Network.u(:,3),'b-','LineWidth',2.5);
quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.g(:,1),Network.g(:,2),Network.g(:,3),'b--','LineWidth',2.5);

hold off

xlabel('x'); ylabel('y'); zlabel('z');
axis equal
return
%%

fun = @(t,y) derivatives(t,y,Network,constants);

fun = @(t,y) derivatives(t,y, Mesh, Network, mesh_node_parameters, index_mapping, matrix_props, assembly_input);

y0 = [ Network.l; reshape(Network.nodes',Network.n_nodes*3,1) ];
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
l = y(:,1:Network.n_links);
for i = 1:size(y,1)
% x(:,:,i) = reshape( y(i,Network.n_links + 1 : Network.n_links + Network.n_nodes*3) , 3,[])';
x = reshape( y(i,Network.n_links + 1 : Network.n_links + Network.n_nodes*3) , 3,[])';

figure(123);
plot3(x(:,1),x(:,2),x(:,3),'ko','MarkerFaceColor','k');
% text(Network.nodes(:,1)+0.2,Network.nodes(:,2),Network.nodes(:,3),cellfun(@num2str,num2cell(1:Network.n_nodes),'UniformOutput',false),'FontSize',12);
hold on
plot3([x(Network.links(:,1),1) x(Network.links(:,2),1)]',...
    [x(Network.links(:,1),2) x(Network.links(:,2),2)]',...
    [x(Network.links(:,1),3) x(Network.links(:,2),3)]',...
    'r--','LineWidth',2);

hold off

xlabel('x'); ylabel('y'); zlabel('z');
axis equal
xlim([-25 25]);  ylim([-25 25]);  zlim([-25 25]); 
title(num2str(t(i)));
% pause(0.1)
drawnow

end