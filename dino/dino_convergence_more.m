clear sols


% uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_medium_Body-Transverse-Coplanar_Hairs_interpolation.fig',1)
% uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_medium_Body-Transverse-Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig',1)
% uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_medium_Transverse-Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig',1)
% uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_medium_Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig',1)
% uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_medium_Normal_Top_Hairs_interpolation.fig',1)
% uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_medium_Normal_Bottom_Hairs_interpolation.fig',1)

%   uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_all_refined_Normal_Bottom_Hairs_interpolation.fig',1)
%  uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_medium_Body-Transverse-Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig',1)
% 

% for Transverse, Coplanar, Normal Top, Normal Bottom (each alone),
% difference from refined < 5%
% refined case is so refined that can only run each thing by itself, so no
% way to test error of entire geom
uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_base2_Transverse_interpolation.fig',1)



chil = get(gcf,'Children');
        for j = 1:6
            sols(1,j) = chil(j).Children.YData;
        end
        
        close
        
%         uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_all_refined_Body-Transverse-Coplanar_Hairs_interpolation.fig',1)
%          uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_all_refined_Body-Transverse-Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig',1)
%         uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_all_refined_Transverse-Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig',1)
%         uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_all_refined_Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig',1)
%          uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_all_refined_Normal_Top_Hairs_interpolation.fig',1)
%          uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_all_refined_more_Normal_Bottom_Hairs_interpolation.fig',1)
%             uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_body_all_refined_Body-Transverse-Normal_Top_Hairs-Normal_Bottom_Hairs_interpolation.fig',1)
     uiopen('C:\Users\rudi\Desktop\RD\dino_dumps_uberhairs\meshes_refined_Transverse_interpolation.fig',1)
%     
        
        chil = get(gcf,'Children');
        for j = 1:6
            sols(2,j) = chil(j).Children.YData;
        end
        
        close
        
        % order of subplots is reversed here....
        U_mag = sqrt(sum(sols(:,4:6).^2,2));
        Omega_mag = sqrt(sum(sols(:,1:3).^2,2));
        
        diff(U_mag) / U_mag(end) * 100
        diff(Omega_mag) / Omega_mag(end) * 100
        
        