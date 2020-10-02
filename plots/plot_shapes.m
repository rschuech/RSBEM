folder = 'E:\Hull\all_meshes\';
outfolder = 'E:\Hull\shapes2\';

bodies = dir([folder,'curved_rod*metadata.mat']);
bodies = {bodies.name};

bad = [];
for i = 1:length(bodies)
    load([folder,bodies{i}]);
    if isfield(Metadata.mesh,'meshing_succeeded') && ~Metadata.mesh.meshing_succeeded
        bad(end+1) = i;
    end
end

bodies(bad) = [];


tails = dir([folder,'tail*metadata.mat']);
tails = {tails.name};

bad = [];
for i = 1:length(tails)
    load([folder,tails{i}]);
    if isfield(Metadata.mesh,'meshing_succeeded') && ~Metadata.mesh.meshing_succeeded
        bad(end+1) = i;
    end
end

tails(bad) = [];


clear AR1 AR2
for i = 1:length(bodies)
    AR1(i) = str2double(bodies{i}(strfind(bodies{i},'AR1_')+4:strfind(bodies{i},'_AR2_')-1));
    AR2(i) = str2double(bodies{i}(strfind(bodies{i},'AR2_')+4:strfind(bodies{i},'_metadata')-1));
end

clear amp lambda nlambda
for i = 1:length(tails)
    amp(i) = str2double(tails{i}(strfind(tails{i},'amp_')+4:strfind(tails{i},'_lambda_')-1));
    lambda(i) = str2double(tails{i}(strfind(tails{i},'lambda_')+7:strfind(tails{i},'_nlambda_')-1));
    nlambda(i) = str2double(tails{i}(strfind(tails{i},'nlambda_')+8:strfind(tails{i},'_metadata')-1));
end


%%
% AR1s = [1  1.5 2 4 6 8 10 12];
% AR2s = [0  0.2 0.4 0.6 0.8 0.9 0.95];
% amps = [0.3015  0.402 0.5025 0.7035  ];
% lambdas = [ 2.1774   2.9032 3.6291  5.0807 ];
% nlambdas = [1.1175   1.49  1.8625];


AR1s = X_Y_unq(:,1);  AR2s = X_Y_unq(:,2);
amps = best_amps;  lambdas = best_lambdas;  nlambdas = best_nlambdas;

runtype = 'plots';
%runtype = 'limits';
%limits.x = [0 0];  limits.y = [0 0];  limits.z = [0 0];
c = 0;
% for i = 1:length(AR1s)
%     for j = 1:length(AR2s)
%         for k = 1:length(amps)
%             for l = 1:length(lambdas)
%                 for m = 1:length(nlambdas)
for i = 1:length(AR1s)
    c = c+1;
    
    %                     c / (length(AR1s)*length(AR2s)*length(amps)*length(lambdas)*length(nlambdas))
    c / length(AR1s)
    
    %                     body_name = ['curved_rod_AR1_',num2str(AR1s(i)),'_AR2_',num2str(AR2s(j))];
    %                     tail_name = ['tail_radius_0.031018_amp_',num2str(amps(k)),'_lambda_',num2str(lambdas(l)),'_nlambda_',num2str(nlambdas(m))];
    
    
    body_name = ['curved_rod_AR1_',num2str(AR1s(i)),'_AR2_',num2str(AR2s(i))];
    tail_name = ['tail_radius_0.031018_amp_',num2str(amps(i)),'_lambda_',num2str(lambdas(i)),'_nlambda_',num2str(nlambdas(i))];%
    try
        [body,body_meta] = load_mesh([folder,body_name,'.dat']);
        [tail,tail_meta] = load_mesh([folder,tail_name,'.dat']);
        %                     catch
        %                         continue
        %                     end
        
        if isempty(body) || isempty(tail)
            continue
        end
        
        body.name = 'body';  tail.name = 'tail';
        body.orientation = [1 0 0; 0 1 0; 0 0 1;]';  %initial orientation vectors are along x, y, z axes
        tail.orientation = [1 0 0; 0 1 0; 0 0 1;]';
        
        body.refpoints = [0 0 0]';
        tail.refpoints = [0 0 0]';  %arbitrary point near tail volume centroid (for orientation vector plotting)
        clear Mesh Metadata
        Mesh(1) = body;  Metadata(1) = body_meta;
        Mesh(2) = tail;  Metadata(2) = tail_meta;
        
        [Mesh] = global_inds(Mesh);
        [Mesh] = renumber_Mesh(Mesh);
        % create smaller input struct with just necessary parameters to reduce mex
        % pain
        temp_input.performance.nthreads = 8;
        temp_input.accuracy.integration_tol.area.abstol = 1E10 ;
        temp_input.accuracy.integration_tol.area.reltol = 1E10 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        temp_input.accuracy.integration_tol.area.maxevals = 1E10;
        temp_input.accuracy.integration_tol.centroid.abstol = 1E10 ;
        temp_input.accuracy.integration_tol.centroid.reltol = 1E10 ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        temp_input.accuracy.integration_tol.centroid.maxevals = 1E10;
        temp_input.accuracy.integration_tol.volume.abstol = 1E10 ;
        temp_input.accuracy.integration_tol.volume.reltol = 1E10  ; %even this is probably way overkill - just one adaptive iteration seems to be extremely accurate
        temp_input.accuracy.integration_tol.volume.maxevals = 1E10;
        
        [Mesh] = store_mesh_constants_wrapper(Mesh, temp_input);
        
        
        tail_radius = 0.031018;
        if AR1s(i) == 1 %sphere degenerate case
            Mesh(2) = shiftMesh(Mesh(2),[-2*Metadata(2).geom.pipeRadius - Metadata(1).geom.radius 0 0]);
        elseif AR2s(i) == 0  %not a sphere, but a straight rod degenerate case
            Mesh(2) = shiftMesh(Mesh(2),[-2*Metadata(2).geom.pipeRadius - (Metadata(1).geom.radius + Metadata(1).geom.height/2) 0 0]);
        else  %actually a general curved rod
            Mesh(2) = shiftMesh(Mesh(2),[-2*Metadata(2).geom.pipeRadius - (Metadata(1).geom.radius) 0 0]);
        end
        
        switch runtype
            case 'limits'
                verts = vertcat(Mesh.verts);
                limits.x = [min(limits.x(1),min(verts(:,1))), max(limits.x(2),max(verts(:,1)))];
                limits.y = [min(limits.y(1),min(verts(:,2))), max(limits.y(2),max(verts(:,2)))];
                limits.z = [min(limits.z(1),min(verts(:,3))), max(limits.z(2),max(verts(:,3)))];
            case 'plots'
                
                figure(123)
                clf
                [s,e] = plot_mesh(Mesh);
                xlim(limits.x);  ylim(limits.y);  zlim(limits.z);
                axis off
                light
                set(gca,'view',[-17 40]);
                set(gcf,'position',[ 680         211        1011         767]);
                drawnow
                print([outfolder,body_name,'_',tail_name,'.png'],'-dpng','-r450');
                %  print([outfolder,body_name,'_',tail_name,'.eps'],'-depsc');
        end
        
        
    end
end



