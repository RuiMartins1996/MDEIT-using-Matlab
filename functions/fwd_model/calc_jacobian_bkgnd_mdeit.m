function jacobian_bkgnd = calc_jacobian_bkgnd_mdeit(img)
    % TODO
    
    jacobian_bkgnd = mean(img.elem_data)*ones(numel(img.elem_data),1);
end