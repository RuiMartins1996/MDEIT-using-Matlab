function r = calc_residual_mdeit(img,x,data,recon_mode,select_sensor_axis)

valid_modes = {'mdeit1', 'mdeit3'};

if ~ismember(recon_mode, valid_modes)
    error('my_inv_solve: invalid recon_mode "%s". Must be''mdeit1'', or ''mdeit3''.', recon_mode);
end

% EIT solve
img.elem_data = x;
u = fwd_solve(img);

if strcmp(recon_mode,'mdeit1')
    s = fwd_solve_mdeit(img,u);
    
    switch select_sensor_axis
        case 1
            r = s.Bx(:)-data(:);
        case 2
            r = s.By(:)-data(:);
        case 3
            r = s.Bz(:)-data(:);
        otherwise
            error('select_sensor_axis must be 1,2 or 3');
    end
    
    return
else
    s = fwd_solve_mdeit(img,u);
    
    r = [s.Bx(:);s.By(:);s.Bz(:)]-data(:);
    return
end

end

