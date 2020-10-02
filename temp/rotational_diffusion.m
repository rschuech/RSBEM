function [    cos_theta ] = rotational_diffusion(Nt,Np,dt,principal_axes, swimming_axis0,D)

cos_theta = NaN(Nt,Np);

parfor j = 1:Np % particle
    current_axes = principal_axes;
    swimming_axis = swimming_axis0;
    cos_theta_temp = NaN(Nt,1); cos_theta_temp(1) = 1;
    for i = 2:Nt % time
        %     i/size(dx,1)
        dx = randn(3,1) .* sqrt(2*D*dt);
        
        order = randperm(length(D));
        % order = [1 2 3];
        for k = 1:length(order)
            %         rotmat = rotate_arbitrary_vector( current_axes(:,order(1)), dx(order(1))) * rotate_arbitrary_vector( current_axes(:,order(2)), dx(order(2))) * rotate_arbitrary_vector( current_axes(:,order(3)), dx(order(3)));
            
            rotmat = rotate_arbitrary_vector( current_axes(:,order(k)), dx(order(k)));
            
            
            current_axes = rotmat * current_axes;   current_axes = current_axes ./ repmat( sqrt(sum(current_axes.^2)) , 3,1);
            swimming_axis = rotmat * swimming_axis;   swimming_axis = swimming_axis ./ sqrt(sum(swimming_axis.^2));
            
        end
        
        cos_theta_temp(i) =   sum( swimming_axis .* swimming_axis0 );
    end
    
    cos_theta(:,j) = cos_theta_temp;
    
end