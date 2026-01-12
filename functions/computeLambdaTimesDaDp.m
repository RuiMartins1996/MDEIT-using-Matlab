function ltdAdsigma = computeLambdaTimesDaDp(img,lambda)

lambda_times_Gx = img.fwd_model.G.Gx*lambda;
lambda_times_Gy = img.fwd_model.G.Gy*lambda;
lambda_times_Gz = img.fwd_model.G.Gz*lambda;

wx = lambda_times_Gx .* img.fwd_model.elem_volume;          % K x 1
wy = lambda_times_Gy .* img.fwd_model.elem_volume;
wz = lambda_times_Gz .* img.fwd_model.elem_volume;

% Scale each row k of G^alpha by w_alpha(k):
Gxw = img.fwd_model.G.Gx .* wx;            % broadcasting: KxN .* Kx1 -> KxN
Gyw = img.fwd_model.G.Gy .* wy;
Gzw = img.fwd_model.G.Gz .* wz;

% Now transpose rows->columns to get N x K matrices and sum:
ltdAdsigma = Gxw' + Gyw' + Gzw';    % N x K

end