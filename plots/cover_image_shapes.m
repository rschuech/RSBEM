
files = dir(['C:\Users\rudi\Desktop\RD\cover image dumps\*torque_dump*.mat']);
files = {files.name};

lims = [Inf -Inf; Inf -Inf; Inf -Inf];
for ff = 1:length(files)
    load(['C:\Users\rudi\Desktop\RD\cover image dumps\', files{ff}]);
    
i = 1;

f = reshape( Solutions.f{i} , 3 , [] )';

F{1} = f(Mesh(1).indices.glob.vert,:);
F{2} = f(Mesh(2).indices.glob.vert,:);

F_mag{1} = sqrt(sum(F{1}.^2,2));
F_mag{2} = sqrt(sum(F{2}.^2,2));

figure(400+ff); clf;
[s,e] = plot_mesh(Mesh,[4 2],F_mag);
axis off

% dirs = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1;];
% for d = 1:size(dirs,1)
%     ligs(d) = light;
%     ligs(d).Position = dirs(d,:);  % [-1 0 1];
% end
% material dull
% set(s,'ambientstrength',0.7);

% pause
temp = [xlim; ylim; zlim];
lims = [min(temp(:,1),lims(:,1))  max(temp(:,2),lims(:,2))];

end

for ff = 1:length(files)
figure(400+ff);
xlim(lims(1,:)); ylim(lims(2,:)); zlim(lims(3,:));

saveas(gcf,['C:\Users\rudi\Desktop\RD\cover image dumps\',files{ff}(1:end-24),'.fig']);
end
