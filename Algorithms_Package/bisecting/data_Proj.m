%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Calculates the projections of the data set V onto the lines chosen
% and returns the projected points as vectors in the cell named
% proj_Cell. Pca_trigger decides whether pca projections of the data
% will be used, axis_trigger whether axis projections will be used,
% randproj_trigger whether random projections will be used and
% randproj_number how many of them.
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[proj_Cell] = data_Proj(V,pca_trigger,axis_trigger,randproj_trigger,randproj_number)

    proj_Cell={};
    chosen = 0;
    
    % pca projections
    if pca_trigger == 1
        % Instead of the first dominant pca component, calculate
        % the pca components that their eigenvalue sum is less
        % than 90 percent of the total eigenvalue sum.
        [pc,latent,~] = pcacov(cov(V));
        lat_sum = sum(latent);
        chosen = 1;
        current_sum = 0;
        
        for i = 1:length(latent)
            current_sum = current_sum + latent(i);
            if current_sum >= 0.9*lat_sum
                break;
            end
            chosen = chosen +1;
        end
 
        for c = 1 : chosen
            pdip_projected = pc(:,c)'*V';
            proj_Cell{c} = pdip_projected;
        end
    end
    
    col_Of_Data = size(V,2);
    
    % axis projections
    if axis_trigger == 1
        proj_Vector = zeros(col_Of_Data,1);

        for i = 1:col_Of_Data
            proj_Vector(i)=1;
            pdip_axis_projected = proj_Vector'*V';
            proj_Cell{i+chosen} = pdip_axis_projected;
            proj_Vector(i)=0;
        end
        axis_trigger = col_Of_Data;
    end

   % random projections
   if randproj_trigger == 1
        for j=1:randproj_number 
            rand_proj_Vector = rand(col_Of_Data,1);
            rand_proj_Vector = rand_proj_Vector./ sqrt(rand_proj_Vector' * rand_proj_Vector);
            pdip_rand_projected = rand_proj_Vector' * V';
            proj_Cell{axis_trigger+chosen+j} = pdip_rand_projected;
       end
   end
   
end