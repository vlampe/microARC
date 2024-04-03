% helper function to interpolate mod
% mod must have size [depth doy traj], depth is first dim
% interpolates from model depths depthsMOD to target depths depthsTARGET
% increases size(mod) from [n doy traj] to [m doy traj]
% fixed :)

function out = interpolate_mod(mod, depthsMOD, depthsTARGET)
    % depthsMOD = vector that holds depth levels of mod (summer.FixedParams.z) [n x 1 double]
    % depthsTARGET = vector that holds target depths (1:sum(zwidth)) [m x 1  double]
        
    modINT = NaN(length(depthsTARGET), size(mod, 2), size(mod, 3));  
    for i = 1:size(mod,2)
        for j = 1:size(mod,3)
            modINT(:, i, j) = interp1(depthsMOD, mod(:,i,j), depthsTARGET); 
        end
    end
    out = modINT;
end