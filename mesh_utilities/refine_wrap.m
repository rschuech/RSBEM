function [Surface_vis] = refine_wrap(Mesh, imesh,  n_refines, Surface_vis) %since Mesh is an input argument and Surface_vis an output, they are local variables
coder.extrinsic('roundn');

nplaces = 9;  %for rounding away roundoff errors

Surface_vis.nodes = Mesh.nodes;
Surface_vis.elements = NaN(4*Mesh.n_elements,3);
Surface_vis.locals = NaN(4*Mesh.n_elements,6);
Surface_vis.xi_eta = NaN(Mesh.n_nodes,2);
Surface_vis.xi_eta_elements = NaN(Mesh.n_nodes,1);
Surface_vis.n_nodes = NaN;
Surface_vis.n_elements = NaN;


for elements_i = 1:Mesh.n_elements
    
    inds = Mesh.elements(elements_i,:);
    new1 = inds([1 4 6]);
    new2 = inds([4 2 5]);
    new3 = inds([6 5 3]);
    new4 = inds([4 5 6]);
    
    Surface_vis.xi_eta(inds,:) = [ 0 0; 1 0; 0 1; 1/2 0; 1/2 1/2; 0 1/2; ];
    Surface_vis.xi_eta_elements(inds) = elements_i;
    
    Surface_vis.elements(1 + (elements_i-1)*4 :1 + (elements_i-1)*4 + 3,:) = [new1; new2; new3; new4 ];
    Surface_vis.locals(1 + (elements_i-1)*4 :1 + (elements_i-1)*4 + 3,:) = [ [0 0 1/2 1/2 1; 1/2 0 1 1/2 1; 0 1/2 1/2 1 1; 0 0 1/2 1/2 2] repmat(elements_i,4,1) ];
    %[min_xi min_eta max_xi max_eta triangle_type orig_elements_number]
end





for n = 1:(n_refines(imesh)-1)
    Surface_vis.nodes(end+1:end+3*size(Surface_vis.elements,1),:) = NaN;
    if isfield(Surface_vis,'surfun')
        Surface_vis.surfun(end+1:end+3*size(Surface_vis.elements,1),:) = NaN;
    end
    nodeind = find(isnan(Surface_vis.nodes(:,1)),1,'first');
    nodeind = nodeind(1);  %Coder is dumb
    Surface_vis.xi_eta(end+1:end+3*size(Surface_vis.elements,1),:) = NaN;
    Surface_vis.xi_eta_elements(end+1:end+3*size(Surface_vis.elements,1),:) = NaN;
    
    newelements = NaN(4*size(Surface_vis.elements,1),3);
    newlocals = NaN(4*size(Surface_vis.elements,1),6);
    newc = 1;
    
    newinds = NaN(1,3);
    
    
    for elements_i = 1:size(Surface_vis.elements,1)
        
        orig_elements = Surface_vis.locals(elements_i,end); %original elementsent index for current triangle
        corners = Surface_vis.elements(elements_i,:);  %corner node inds for current triangle
        
        subnodes = Mesh.nodes(Mesh.elements(orig_elements,:),:); %orig quad elementsent nodes containing this triangle
        
        if isfield(Surface_vis,'surfun')
            subsurfun = Surface_vis.surfun(Mesh.elements(orig_elements,:),:);
        end
        
        shape_parameters = Mesh.shape_parameters(orig_elements,:);  %[alpha beta gamma]
        
        coords = Surface_vis.locals(elements_i,1:4);
        
        x1 = coords(1);  x2 = coords(3);  y1 = coords(2);  y2 = coords(4);
        xavg = (x1+x2)/2;  yavg = (y1+y2)/2;
        
        
        
        switch Surface_vis.locals(elements_i,5)
            case 1 %normal triangle
                
                %1
                newlocals(newc,:) = [x1 y1 xavg yavg 1 orig_elements];
                xi = xavg;  eta = y1;
                Surface_vis.xi_eta(nodeind,:) = [xi eta];
                Surface_vis.xi_eta_elements(nodeind) = orig_elements;
                [x_interp2, phi] = T6interp(subnodes,xi,eta,shape_parameters);
                if isfield(Surface_vis,'surfun')
                    surfun_interp2 = subsurfun' * phi;
                end
                xi = x1;  eta = yavg;
                Surface_vis.xi_eta(nodeind+1,:) = [xi eta];
                Surface_vis.xi_eta_elements(nodeind+1) = orig_elements;
                [x_interp3, phi] = T6interp(subnodes,xi,eta,shape_parameters);
                if isfield(Surface_vis,'surfun')
                    surfun_interp3 = subsurfun' * phi;
                end
                Surface_vis.nodes(nodeind:nodeind+1,:) = [x_interp2'; x_interp3']; %2 new mid nodes but 1st corner is already there
                if isfield(Surface_vis,'surfun')
                    Surface_vis.surfun(nodeind:nodeind+1,:) = [surfun_interp2'; surfun_interp3'];
                end
                newinds([1 3]) = nodeind:nodeind+1; %size(Surface_vis.nodes,1)-1:size(Surface_vis.nodes,1);  %goes with 2 new mid nodes
                newelements(newc,:) = [corners(1) newinds(1) newinds(3)];
                
                %2
                newlocals(newc+1,:) = [xavg y1 x2 yavg 1 orig_elements];
                
                xi = xavg;  eta = yavg;
                [x_interp3, phi] = T6interp(subnodes,xi,eta,shape_parameters);
                if isfield(Surface_vis,'surfun')
                    surfun_interp3 = subsurfun' * phi;
                end
                Surface_vis.xi_eta(nodeind+2,:) = [xi eta];
                Surface_vis.xi_eta_elements(nodeind+2) = orig_elements;
                Surface_vis.nodes(nodeind+2,:) = [x_interp3']; %1 new mid node
                if isfield(Surface_vis,'surfun')
                    Surface_vis.surfun(nodeind+2,:) = [ surfun_interp3'];
                end
                newinds(2) = nodeind+2;  %last new mid node
                
                newelements(newc+1,:) = [ newinds(1) corners(2) newinds(2)];
                
                %3
                newlocals(newc+2,:) = [x1 yavg xavg y2 1 orig_elements];
                
                newelements(newc+2,:) = [ newinds(3)  newinds(2) corners(3)];
                
                %4
                newlocals(newc+3,:) = [x1 y1 xavg yavg 2 orig_elements];
                
                newelements(newc+3,:) = [ newinds(1)  newinds(2) newinds(3)];
                
            case 2 %flipped triangle
                
                
                %1
                newlocals(newc,:) = [xavg y1 x2 yavg 2 orig_elements];
                xi = x2;  eta = yavg;
                [x_interp2, phi] = T6interp(subnodes,xi,eta,shape_parameters);
                if isfield(Surface_vis,'surfun')
                    surfun_interp2 = subsurfun' * phi;
                end
                Surface_vis.xi_eta(nodeind,:) = [xi eta];
                Surface_vis.xi_eta_elements(nodeind) = orig_elements;
                xi = xavg;  eta = yavg;
                [x_interp3, phi] = T6interp(subnodes,xi,eta,shape_parameters);
                if isfield(Surface_vis,'surfun')
                    surfun_interp3 = subsurfun' * phi;
                end
                Surface_vis.xi_eta(nodeind+1,:) = [xi eta];
                Surface_vis.xi_eta_elements(nodeind+1) = orig_elements;
                
                Surface_vis.nodes(nodeind:nodeind+1,:) = [x_interp2' ; x_interp3']; %2 new mid nodes but 1st corner is already there
                if isfield(Surface_vis,'surfun')
                    Surface_vis.surfun(nodeind:nodeind+1,:) = [surfun_interp2'; surfun_interp3'];
                end
                newinds([2 1]) = nodeind:nodeind+1;  %goes with 2 new mid nodes
                
                newelements(newc,:) = [corners(1) newinds(2) newinds(1)];
                
                %2
                newlocals(newc+1,:) = [xavg yavg x2 y2 2 orig_elements];
                xi = xavg;  eta = y2;
                [x_interp3] = T6interp(subnodes,xi,eta,shape_parameters);
                if isfield(Surface_vis,'surfun')
                    surfun_interp3 = subsurfun' * phi;
                end
                Surface_vis.xi_eta(nodeind+2,:) = [xi eta];
                Surface_vis.xi_eta_elements(nodeind+2) = orig_elements;
                Surface_vis.nodes(nodeind+2,:) = [x_interp3']; %1 new mid node
                if isfield(Surface_vis,'surfun')
                    Surface_vis.surfun(nodeind+2,:) = [ surfun_interp3'];
                end
                newinds(3) = nodeind+2;  %last new mid node
                
                newelements(newc+1,:) = [ newinds(2) corners(2) newinds(3)];
                
                %3
                newlocals(newc+2,:) = [x1 yavg xavg y2 2 orig_elements];
                
                newelements(newc+2,:) = [ newinds(1)  newinds(3) corners(3)];
                
                %4
                newlocals(newc+3,:) = [xavg yavg x2 y2 1 orig_elements];
                
                newelements(newc+3,:) = [ newinds(1)  newinds(2) newinds(3)];
        end
        nodeind = nodeind + 3;
        newc = newc + 4;
    end
    
    
    Surface_vis.elements = newelements;
    Surface_vis.locals = newlocals;
    
    
end

%remove replicated nodeices (better to not replicate them in the first
%place, but this is easier to code....)
ic = NaN(size(Surface_vis.nodes,1),1);
ia = ic;
[un,ia,ic] = unique(roundn(Surface_vis.nodes,-nplaces),'rows');  %hard coded tolerance here, seems to work so far...
% ia = ia(:,1);
%un = A(ia,:) and A = un(ic,:)
% ic = ic(:,1); %%Coder is dumb

temp = Surface_vis.elements;
for i = 1:size(Surface_vis.elements,1)
    %Surface_vis.elements(i,:);
    temp(i,:) = ic(Surface_vis.elements(i,:));
end

Surface_vis.nodes = un;
if isfield(Surface_vis,'surfun')
    Surface_vis.surfun = Surface_vis.surfun(ia,:);
end
Surface_vis.xi_eta = Surface_vis.xi_eta(ia,:);
Surface_vis.xi_eta_elements = Surface_vis.xi_eta_elements(ia);
Surface_vis.elements = temp;

Surface_vis.n_nodes = size(Surface_vis.nodes,1);
Surface_vis.n_elements = size(Surface_vis.elements,1);

end