m = 1;

if exist('fn')
for i = 1:length(fn)
    delete(fn{i});
end
end
if exist('ft')
for i = 1:length(ft)
    delete(ft{i});
end
end
if exist('fq')
for i = 1:length(fq)
    delete(fq{i});
end
end

fn = [];
ft = [];
fq = [];


for i = 1:Mesh(m).n_vert
    
    normals =  Mesh(m).normals_vert{i};
    if isempty(normals)
        error('shit')
        
    end
    
    inds = combnk(1:size(normals,1),2);
    
    angles = [];
    for ii = 1:size(inds,1)
        
        angles(ii) =   acos( dot( normals(inds(ii,1),:) , normals(inds(ii,2),:) ) ) * 180/pi;
        
    end
    
    if  true || Mesh(m).vert_type(i) == 2 %   true || max(angles) > 10
%              try, delete(fq), end;
        %      try, delete(fq2), end
%         clear fqfq2
        for ii = 1:size(normals,1)
            hold on;
%                      fq{end+1} = quiver3(Mesh(m).verts(i,1),Mesh(m).verts(i,2),Mesh(m).verts(i,3),normals(ii,1),normals(ii,2),normals(ii,3)); 
%                      fq{end}.Color = 'b'; fq{end}.LineWidth = 2; fq{end}.MaxHeadSize = 10;
        end
        
       fn{end+1} = quiver3(Mesh(m).verts(i,1),Mesh(m).verts(i,2),Mesh(m).verts(i,3),Mesh(m).normals_avg(i,1),Mesh(m).normals_avg(i,2),Mesh(m).normals_avg(i,3));  
         fn{end}.Color = 'k'; fn{end}.LineWidth = 2; fn{end}.MaxHeadSize = 10;
      
%          tangents =  eye(3) - Mesh(m).normals_avg(i,:)' * Mesh(m).normals_avg(i,:);
%        tangents = tangents ./ sqrt(sum(tangents.^2,2));
%        
%        dots = [dot(tangents(3,:),tangents(2,:)); dot(tangents(1,:),tangents(3,:)); dot(tangents(2,:),tangents(1,:))];
%        [~,ind] = min(abs(dots));
%        tangents(ind,:) = [];
       
       for j = 1:2
        ft{end+1} = quiver3(Mesh(m).verts(i,1),Mesh(m).verts(i,2),Mesh(m).verts(i,3),Mesh(m).tangents_avg(i,1,j),Mesh(m).tangents_avg(i,2,j),Mesh(m).tangents_avg(i,3,j)); 
       ft{end}.Color = 'b'; ft{end}.LineWidth = 2; ft{end}.MaxHeadSize = 10;
       end
       
%         pause
    end
    
end

%  if ~isempty(angles)
%      angles
%  end
