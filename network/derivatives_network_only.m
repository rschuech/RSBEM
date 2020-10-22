function [dydt] = derivatives(t,y,Stokeslets,constants)

l = y(1:Stokeslets.n_links);
nodes = reshape( y(Stokeslets.n_links + 1 : Stokeslets.n_links + Stokeslets.n_nodes*3) , 3,[])';

d = nodes(Stokeslets.links(:,2),:)  - nodes(Stokeslets.links(:,1),:) ;
r = sqrt(sum( d.^2 , 2 ) );

f_s = Stokeslets.l_0.^2 .* Stokeslets.E .* (r./l - 1) .* d ./ r;


g = NaN(Stokeslets.n_nodes,3); % net Maxwell element viscoelastic force on each network node due to all attached links
u = NaN(Stokeslets.n_nodes,3);

parfor i = 1:Stokeslets.n_nodes
    g(i,:) = sum(   Stokeslets.link_members{i}(2,:)' .*  f_s( Stokeslets.link_members{i}(1,:) , :) , 1);
    
end

parfor i = 1:Stokeslets.n_nodes
    temp = 0;
    for j = 1:Stokeslets.n_nodes
        temp = temp + g(j,:) *  calcS(nodes(j,:)',nodes(i,:)',constants.epsilon^2);
    end
    u(i,:) = 1/constants.mu * 1/constants.alpha * 1/2 * temp;
end


dldt =  Stokeslets.E .* Stokeslets.l_0 ./ Stokeslets.eta .* ( r ./ l - 1);



dydt = [ dldt; reshape(u',Stokeslets.n_nodes*3,1) ];
