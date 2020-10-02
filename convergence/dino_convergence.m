submeshes = {'Body','Transverse','Coplanar_Hairs','Normal_Top_Hairs','Normal_Bottom_Hairs'};
inds = [2 4 5 6];
inds = 1;

fines = {'coarse','medium-coarse','medium','fine','finer'};
fines2 = {'coarse','medium_coarse','medium','fine','finer'};

fines = {'coarse','medium','finer'};
fines2 = {'coarse','medium','finer'};

folder = 'C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\';
clear sols
%   L =  [0.8113  0.5307  0.3603  0.1701  0.0869];

for i = 1%:length(submeshes)
    for ii = 1:length(fines)
        %file = ['meshes_',fines{ii},'_partial_turns2_Body-',submeshes{i},'_interpolation.fig'];
%         file = ['meshes_',fines{ii},'_partial_turns2_',submeshes{i},'_interpolation.fig'];
          file = ['meshes_body_',fines{ii},'_Body-Transverse-Coplanar_Hairs-Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig'];
        
        open([folder,file]);
        
        
        chil = get(gcf,'Children');
        for j = 1:6
            sols.(submeshes{i})(ii,j) = chil(j).Children.YData;
        end
        close
        
    end
end



  %%
 
  for j = 1%:length(submeshes)
      clear L
      for jj = 1:length(fines2)
      L(jj) = mean(sqrt(Meshes.(fines2{jj})(inds(j)).area));
      end
      figure(j+100)
   
      temp = repmat( sols.(submeshes{j})(end,:)  , length(fines), 1 );
      diffs = abs( ( sols.(submeshes{j}) - temp )  ./  temp)  * 100;
     
%       plot(L(1:end-1),diffs(1:end-1,1),'or-','linewidth',2);
%       hold on
%       plot(L(1:end-1),diffs(1:end-1,2),'or--','linewidth',2);
%       plot(L(1:end-1),diffs(1:end-1,3),'or:','linewidth',2);
%          plot(L(1:end-1),diffs(1:end-1,4),'ob-','linewidth',2);
%       plot(L(1:end-1),diffs(1:end-1,5),'ob--','linewidth',2);
%       plot(L(1:end-1),diffs(1:end-1,6),'ob:','linewidth',2);
      
      U_mag = sqrt(sum(sols.(submeshes{j})(:,1:3).^2 , 2));
      temp = repmat( U_mag(end,:)  , length(fines), 1 );
      diffs = abs( (U_mag - temp )  ./  temp)  * 100;
       plot(L(1:end-1),diffs(1:end-1,1),'or-','linewidth',3,'markerfacecolor','r');  hold on;
      
       Omega_mag = sqrt(sum(sols.(submeshes{j})(:,4:6).^2 , 2));
       temp = repmat( Omega_mag(end,:)  , length(fines), 1 );
      diffs = abs( (Omega_mag - temp )  ./  temp)  * 100;
        plot(L(1:end-1),diffs(1:end-1,1),'ob-','linewidth',3,'markerfacecolor','b');
      
      hold off
%       plot(L(1:end-1),diffs(1:end-1,:),'o-'); 
grid on
  end