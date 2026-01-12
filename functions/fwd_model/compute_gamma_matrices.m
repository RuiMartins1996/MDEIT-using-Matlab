function img = compute_gamma_matrices(img)

num_sensors = numel(img.fwd_model.sensors);
num_nodes   = size(img.fwd_model.nodes, 1);
expectedSize = [num_sensors, num_nodes];

mu_factor = img.fwd_model.mu0/(4*pi);

% Convenience handles
R = img.fwd_model.R;
G = img.fwd_model.G;

% Sigma = sparse(1:length(img.elem_data), 1:length(img.elem_data), img.elem_data);
Sigma = spdiags(img.elem_data(:), 0, length(img.elem_data), length(img.elem_data));

% Check for existing fields and reuse if available
if isfield(img, 'GammaX') && ~isempty(img.GammaX)
    if not(isequal(size(img.GammaX), expectedSize))
        error('Expected size for Gamma matrix is wrong')
    end
else
    GammaX = mu_factor * ( -R.Rz * Sigma * G.Gy +  R.Ry * Sigma * G.Gz );
    img.GammaX = GammaX;
end

if isfield(img, 'GammaY') && ~isempty(img.GammaY)
    if not(isequal(size(img.GammaY), expectedSize))
        error('Expected size for Gamma matrix is wrong')
    end
else
    GammaY = mu_factor * ( -R.Rx * Sigma * G.Gz +  R.Rz * Sigma * G.Gx );
    img.GammaY = GammaY;
end

if isfield(img, 'GammaZ') && ~isempty(img.GammaZ)
    if not(isequal(size(img.GammaZ), expectedSize))
        error('Expected size for Gamma matrix is wrong')
    end
else
    GammaZ = mu_factor * ( -R.Ry * Sigma * G.Gx +  R.Rx * Sigma * G.Gy );
    img.GammaZ = GammaZ;
end

end
