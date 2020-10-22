
m = [8 6 6]; % # nodes in each dimension
x = linspace(-2.25,1,m(1)); % cube 20 um wide
y = linspace(-1,1,m(2));
z = linspace(-1,1,m(3));

% single plane of points around tail
% m = [1 12 12];
% x = linspace(-1.75,-1.75,m(1)); % cube 20 um wide
% y = linspace(-1,1,m(2));
% z = linspace(-1,1,m(3));

% single plane of points around body
m = [1 8 8];
x = linspace(0.65,0.65,m(1)); % cube 20 um wide
y = linspace(-1,1,m(2));
z = linspace(-1,1,m(3));


[X,Y,Z] = ndgrid(x,y,z);

% constants.alpha = 4*pi;  % solid angle for 0D nodes is 4 pi since no fluid is blocked by the geometry


Network.n_nodes = numel(X);
Network.nodes = [X(:) Y(:) Z(:)];

% NaN(3*m^2*(m-1),2) % number of links from Wrobel et al 2016
links = NaN(0,2);
sz = size(X);
for s = 1:Network.n_nodes
    [i,j,k] = ind2sub(sz,s);
    
    stencil = [i-1, j, k;  i+1,j,k;  i,j-1,k;  i,j+1,k;  i,j,k-1;  i,j,k+1;];
    for st = 1:size(stencil,1)
        if any(stencil(st,:) < 1) || any(stencil(st,:) > sz) % stencil point doesn't lie within range of actual grid
            continue
        end
        links(end+1,:) = [sub2ind(sz,stencil(st,1),stencil(st,2),stencil(st,3)) s];
        
    end
end

links = unique(sort(links,2),'rows');


 links = zeros(0,2);

Network.n_links = size(links,1);
Network.links = links;

Network.nodes = [ Network.nodes; Mesh(1).nodes(1:2:end,:); Mesh(2).nodes(1:4:end,:)];
Network.n_nodes = size(Network.nodes,1);

Network.E = repmat(1,Network.n_links,1);
Network.l_0 = sqrt(sum( (Network.nodes(Network.links(:,1),:)  - Network.nodes(Network.links(:,2),:) ).^2 , 2 ) );
% Network.l_0 =  sqrt(sum(  (( Network.nodes(Network.links(:,1),:)  + Network.nodes(Network.links(:,2),:) ) / 2).^2 , 2));
% resting lengths are the initial distances from the center of each link to the origin
Network.l = Network.l_0;
Network.eta = repmat(10,Network.n_links,1);
Network.link_members = cell(Network.n_nodes,1);

for i = 1:Network.n_nodes
    Network.link_members{i,1} = NaN(2,0);
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



%%%%% Compute g, the resultant force at each network node due to all viscoelastic links with other network nodes

d = Network.nodes(Network.links(:,2),:)  - Network.nodes(Network.links(:,1),:) ;
r = sqrt(sum( d.^2 , 2 ) );

f_s = Network.l_0.^2 .* Network.E .* (r./Network.l - 1) .* d ./ r;

g = NaN(Network.n_nodes,3); % net Maxwell element viscoelastic force on each network node due to all attached links

parfor i = 1:Network.n_nodes
    g(i,:) = sum(   Network.link_members{i}(2,:)' .*  f_s( Network.link_members{i}(1,:) , :) , 1);
    
end

Network.g = g;  clear g;




%%

figure(124);
plot3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),'ko','MarkerFaceColor','k');
% text(Network.nodes(:,1)+0.2,Network.nodes(:,2),Network.nodes(:,3),cellfun(@num2str,num2cell(1:Network.n_nodes),'UniformOutput',false),'FontSize',12);
hold on
plot3([Network.nodes(Network.links(:,1),1) Network.nodes(Network.links(:,2),1)]',...
    [Network.nodes(Network.links(:,1),2) Network.nodes(Network.links(:,2),2)]',...
    [Network.nodes(Network.links(:,1),3) Network.nodes(Network.links(:,2),3)]',...
    'r--','LineWidth',2);

% quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.u(:,1),Network.u(:,2),Network.u(:,3),'b-','LineWidth',2.5);
% quiver3(Network.nodes(:,1),Network.nodes(:,2),Network.nodes(:,3),Network.g(:,1),Network.g(:,2),Network.g(:,3),'b--','LineWidth',2.5);

[s,e] = plot_mesh(Mesh,2);
light



hold off

xlabel('x'); ylabel('y'); zlabel('z');
axis equal
return
