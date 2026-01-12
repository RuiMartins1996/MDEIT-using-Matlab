function sigma_std = noise_correction(imgh,imgi,Jh,lambda,noisy_data_generator,num_noise_repetitions)

% Solve Tykhonov regularized normal equations using SVD for noisy inputs

[U,S,V] = svd(Jh,'econ');

n_elem = size(Jh,2);
M = zeros(n_elem,num_noise_repetitions);

for t = 1:num_noise_repetitions
    
    fprintf('Solving for noise realization %i\n',t);
    
    data_noisy = noisy_data_generator(imgh,imgi);
    
    sv = diag(S)+lambda./diag(S);
    sigma_noisy = V*diag(1./sv)*U'*data_noisy;
    
    M(:,t) = sigma_noisy;
end

sigma_std = std(M,[],2);

end

