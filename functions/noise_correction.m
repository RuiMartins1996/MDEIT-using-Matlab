function sigma_std = noise_correction(imgh,imgi,Jh,lambda,noisy_data_generator,num_noise_repetitions,U,S,V)

% imgh and imgi are for data generation, so they should correspond to a
% fine model, not a coarse one.

if nargin <7
    [U,S,V] = svd(Jh,'econ');
end

try % noisy data generator generates all samples once
    data_noisy = noisy_data_generator(imgh,imgi,num_noise_repetitions);
    
    s  = diag(S);
    sv = s + lambda./s;
    M = V * diag(1./sv) * U' * data_noisy;

    sigma_std = std(M,[],2);
catch % noisy data generator generates a sample each time
    

    s  = diag(S);
    sv = s + lambda./s;
    
    M = zeros(size(Jh,2),num_noise_repetitions);

    for i = 1:num_noise_repetitions
        data_noisy = noisy_data_generator(imgh,imgi);
        M(:,i) = V * diag(1./sv) * U' * data_noisy;
    end

    sigma_std = std(M,[],2);
end


end


