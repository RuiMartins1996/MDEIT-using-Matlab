clc; clear all; close all;

rng(1)

% Data of EnvironmentalNoise in: 
%   https://zenodo.org/records/7872660

% Article where the authors collected the data: 
%   https://www.sciencedirect.com/science/article/pii/S1053811923004032
%% Prepare workspace
% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

% Have to add the functions path manually so prepare_workspace runs
parent_folder = fileparts(script_folder);
grandparent_folder =fileparts(parent_folder);
addpath(genpath(fullfile(grandparent_folder,'functions')));

ensure_spm12(); % check if spm12 exists locally in the external folder. If not, prompt the user to download it
ensure_noise_sample_file(); % check if the EnvironmentalNoise folder is in the script's directoty.

model_folder = prepare_workspace(script_folder);

%% Create/fetch model_parameters

SNR_vector = [1 10 30 50];
noise_type = 'opm';

noise_data_path = fullfile(script_folder, ...
    'EnvironmentalNoise\sub-001\ses-001\meg\dsub-001_ses-001_task-noise_run-001_meg.mat');

%IMPORTANT!!!!!!!!!!!!!!!!!!!!

% This does not work out of the box form some reason. Need to add this
% snippet in SPM-12:

% %Rui's change
% %--------------------------------------------------------------------------
% pos = [posOri.Px,posOri.Py,posOri.Pz];
% ori = [posOri.Ox,posOri.Oy,posOri.Oz];
% cl = posOri.name;
% 
% grad= [];
% grad.label = cl;
% grad.coilpos = pos;
% grad.coilori = ori;
% grad.tra = eye(numel(grad.label));
% grad.chanunit = repmat({'T'}, numel(grad.label), 1);
% grad.chantype= 'MEG';
% grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'mm');
% D = sensors(D, 'MEG', grad);
% save(D);
% %--------------------------------------------------------------------------

% In the function spm_opm_create from the spm-12 library in the external folder.
% This runs the %-Place Sensors  in object 
% snippet regardless of the state of the variable "forward"

if ~isfile(noise_data_path)
    S = [];
    S.data = 'EnvironmentalNoise\sub-001\ses-001\meg\sub-001_ses-001_task-noise_run-001_meg.bin';
    S.positions = 'EnvironmentalNoise\sub-001\ses-001\meg\positions.tsv';
    S.path = 'EnvironmentalNoise\analysedData';
    fprintf('Loading: will take around 30s ...\n');
    D = spm_opm_create(S);
    fprintf('Downsampling: will take around 10min ...\n')
    % Downsample
    S = [];
    S.D = D;
    S.fsample_new = 1000;
    D = spm_eeg_downsample(S);
end

%Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);

background_conductivity = 3.28e-1/sigma0;  
anomaly_conductivity = 1e-12/sigma0;

maxsz_reconstruction = 0.05;

num_noise_repetitions = 30;

model_parameters = create_kai_2d_model_parameters(l0, z0, sigma0, I0);

%% Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
stimulation = mk_stim_patterns(16,1,inj,meas,{},current_amplitude);

%% Create model

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl = fmdls{1};
fmdl.stimulation = stimulation;

%% Make images from forward models

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);

% Add plastic cylinder
imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

%% Generate coarse forward model for reconstruction (different mesh than the data)
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction = fmdls{1};

% Don't forget to assignt the same stimulation pattern
fmdl_reconstruction.stimulation = stimulation;

n_elem = size(fmdl_reconstruction.elems,1);

%% Make inverse model
imdl_mdeit_1= eidors_obj('inv_model','my_inverse_model');

% 1-axis imdl
imdl_mdeit_1.recon_mode = 'mdeit1';
imdl_mdeit_1.select_sensor_axis = 3; % z-axis only page 172

imdl_mdeit_1.fwd_model = fmdl_reconstruction; %Use a different forward model for reconstruction
imdl_mdeit_1.jacobian_bkgnd = struct('value',background_conductivity);
imdl_mdeit_1.solver = 'gn';
imdl_mdeit_1.RtR_prior = @(x,Jx) x; %tykhonov
imdl_mdeit_1.recon_type = 'difference';

imdl_mdeit_1.verbose = true; % print debug info


imdl_eit = mk_common_model('a2c2',8); % Will replace most fields
imdl_eit.jacobian_bkgnd = struct('value',background_conductivity); %If I
% use this, the solution blows up
imdl_eit.fwd_model = fmdl_reconstruction; %Use a different forward model for reconstruction
imdl_eit.recon_mode = 'eit';
imdl_eit.reconst_type = 'difference';
imdl_eit.RtR_prior = @prior_tikhonov;


%% Generate data

if strcmp(noise_type,'opm')
    fprintf('Noise type opm still means EIT is white noise!!!!!!')
end

for i = 1:length(SNR_vector)

    fprintf('Running with SNR = %i\n',SNR_vector(i));

    options.data_path = noise_data_path;
    options.noise_structure.type = noise_type;
    options.noise_structure.snr = SNR_vector(i);
    options.noise_structure.epoch_time = 50e-3; %50 milliseconds
    options.noise_structure.sampling_rate = 1000;
    options.noise_structure.B0 = 1.2566370612720e-6*J0*l0;
    options.noise_structure.measurement_protocol = 'default_rui';

    [data_b,data_u,snr,noise_b,noise_u] = generate_data(options,imgh,imgi);

    if abs(snr-SNR_vector(i))>0.01 %in case of opm noise, override snr by the one outputed by generate_data
        SNR_vector (i)= snr;
    end

    % Since this is difference noisy data, it outputs the noisy difference
    % measurements, we have to do the following:

    % Forward solve
    [datah,uh] = fwd_solve_mdeit(imgh);
    Bzh = datah.Bz(:);

    % Correction
    Bzi = Bzh + data_b.dBz;

    [~,ui] = fwd_solve_mdeit(imgi);
    ui.meas = uh.meas +data_u.meas_du;

    %% Find optimal regularization parameter with GCV
    lambda_vector = logspace(-17,3,30);

    [lambda_optimal_mdeit1,optimal_id_mdeit_1,V_mu_mdeit1,~] = generalized_cross_validation(imdl_mdeit_1,data_b.dBz,lambda_vector);
    imdl_mdeit_1.hyperparameter = struct('value',lambda_optimal_mdeit1);


    [lambda_optimal_eit,optimal_id_eit,V_mu_eit,~] = generalized_cross_validation(imdl_eit,data_u.meas_du,lambda_vector);
    % Since EIDORS uses hp^2*RtR, but the generalized_cross_validation function
    % uses hp*RtR, we have to correct for that.
    imdl_eit.hyperparameter = struct('value',sqrt(lambda_optimal_eit));

    %% Perform reconstruction

    % Reconstruct (difference) for 1-MDEIT
    img_output_mdeit_1{i} = inv_solve_mdeit(imdl_mdeit_1,Bzh,Bzi);
    img_output_mdeit_1{i}.hyperparameter = imdl_mdeit_1.hyperparameter;

    % Reconstruct (difference) for EIT
    img_output_eit{i} = inv_solve(imdl_eit,uh,ui);
    img_output_eit{i}.hyperparameter = imdl_eit.hyperparameter;

end

%% Show reconstruction
% 
% show_fem_transparent_edges(imgi)
%     box on
%     title('Ground-Truth','Interpreter','latex')
% 
figure('Position',[200,200,1000,500]);
for i = 1:length(SNR_vector)
    subplot(2,length(SNR_vector),i)
    show_fem_transparent_edges(img_output_eit{i})
    eidors_colourbar(img_output_eit{i})
    box on
    title(sprintf('EIT - SNR = %.f - noise - %s',SNR_vector(i),noise_type),'Interpreter','latex')
    subplot(2,length(SNR_vector),length(SNR_vector)+i)
    show_fem_transparent_edges(img_output_mdeit_1{i})
    eidors_colourbar(img_output_mdeit_1{i})
    box on
    title(sprintf('MDEIT1 - SNR = %.f - noise - %s',SNR_vector(i),noise_type),'Interpreter','latex')
end

%% Noise correction
function data_noisy = noisy_data_generator_mdeit(imgh,imgi,options)
[~,~,~,noise_b,~] =  generate_data(options,imgh,imgi);
data_noisy = noise_b.noise_b_z;
end

function data_noisy = noisy_data_generator_eit(imgh,imgi,options)
[~,~,~,~,noise_u] =  generate_data(options,imgh,imgi);
data_noisy = noise_u;
end

for i = 1:length(SNR_vector)

    options.noise_structure.snr = SNR_vector(i);
    
    func_eit = @(imgh,imgi) noisy_data_generator_eit(imgh,imgi,options);
    func_mdeit = @(imgh,imgi) noisy_data_generator_mdeit(imgh,imgi,options);

    lambda_mdeit = img_output_mdeit_1{i}.hyperparameter.value;
    Jh = img_output_mdeit_1{i}.jacobian;

    sigma_std_mdeit{i} = noise_correction(imgh,imgi,Jh,lambda_mdeit,func_mdeit,num_noise_repetitions);
    
    img_output_mdeit_1_nbc{i} = img_output_mdeit_1{i};
    img_output_mdeit_1_nbc{i}.elem_data = img_output_mdeit_1_nbc{i}.elem_data./sigma_std_mdeit{i};
    
    lambda_eit = img_output_eit{i}.hyperparameter.value^2;
    img_eit_h = mk_image(fmdl_reconstruction,background_conductivity);
    Jh = calc_jacobian(img_eit_h);
    sigma_std_eit{i} = noise_correction(imgh,imgi,Jh,lambda_eit,func_eit,num_noise_repetitions);
    
    img_output_eit_nbc{i} = img_output_eit{i};
    img_output_eit_nbc{i}.elem_data = img_output_eit_nbc{i}.elem_data./sigma_std_eit{i};
end

%% No noise correction
figure('Position',[200,200,1000,500]);

for i = 1:length(SNR_vector)
    subplot(4,length(SNR_vector),i)
    show_fem_transparent_edges(img_output_eit{i})
    eidors_colourbar(img_output_eit{i})
    box on
    title(sprintf('EIT - SNR = %.f - noise - %s',SNR_vector(i),noise_type),'Interpreter','latex')
    subplot(4,length(SNR_vector),1*length(SNR_vector)+i)
    show_fem_transparent_edges(img_output_mdeit_1{i})
    eidors_colourbar(img_output_mdeit_1{i})
    box on
    title(sprintf('MDEIT1 - SNR = %.f - noise - %s',SNR_vector(i),noise_type),'Interpreter','latex')
end

%% With noise correction
for i = 1:length(SNR_vector)
    subplot(4,length(SNR_vector),2*length(SNR_vector)+i)
    show_fem_transparent_edges(img_output_eit_nbc{i})
    eidors_colourbar(img_output_eit_nbc{i})
    box on
    title(sprintf('EIT - SNR = %.f - noise - %s',SNR_vector(i),noise_type),'Interpreter','latex')
    subplot(4,length(SNR_vector),3*length(SNR_vector)+i)
    show_fem_transparent_edges(img_output_mdeit_1_nbc{i})
    eidors_colourbar(img_output_mdeit_1_nbc{i})
    box on
    title(sprintf('MDEIT1 - SNR = %.f - noise - %s',SNR_vector(i),noise_type),'Interpreter','latex')
end


%% Metrics NO NBC

for i = 1:length(SNR_vector)
    
    metrics_eit = compute_greit_metrics_2d(img_output_eit{i},128);
    metrics_mdeit = compute_greit_metrics_2d(img_output_mdeit_1{i},128);
    
    centroid_diff_eit = norm(metrics_eit.centroid(:)-model_parameters.anomaly.position(1:2,1),2);
    centroid_diff_mdeit = norm(metrics_mdeit.centroid(:)-model_parameters.anomaly.position(1:2,1),2);

    fprintf('SNR = %i:\n',SNR_vector(i));
    fprintf('EIT centroid : (%.4f,%.4f) | radius : %.4f | centroid difference:  %.4f\n',metrics_eit.centroid(1),metrics_eit.centroid(2),metrics_eit.radius,centroid_diff_eit);
    fprintf('MDEIT centroid : (%.4f,%.4f) | radius : %.4f | centroid difference:  %.4f\n',metrics_mdeit.centroid(1),metrics_mdeit.centroid(2),metrics_mdeit.radius,centroid_diff_mdeit);
end
%% Metrics with NBC
for i = 1:length(SNR_vector)
    
    metrics_eit = compute_greit_metrics_2d(img_output_eit_nbc{i},128);
    metrics_mdeit = compute_greit_metrics_2d(img_output_mdeit_1_nbc{i},128);
    
    anomaly_true_position = model_parameters.anomaly.position(1:2);
    anomaly_true_radius = model_parameters.anomaly.radius;

    centroid_diff_eit = norm(metrics_eit.centroid(:)-anomaly_true_position(:),2);
    centroid_diff_mdeit = norm(metrics_mdeit.centroid(:)-anomaly_true_position(:),2);

    radius_diff_eit = abs(metrics_eit.radius-anomaly_true_radius);
    radius_diff_mdeit = abs(metrics_mdeit.radius-anomaly_true_radius);

    fprintf('-----------------------------------------\n');
    fprintf('SNR = %2.1f: \n',SNR_vector(i));
    fprintf(['EIT centroid : (%.4f,%.4f) | radius : %.4f \n' ...
        'centroid difference:  %.4f, radius difference: %.4f \n'],...
        metrics_eit.centroid(1),metrics_eit.centroid(2),metrics_eit.radius,centroid_diff_eit,radius_diff_eit);
    fprintf(['MDEIT centroid : (%.4f,%.4f) | radius : %.4f \n' ...
        'centroid difference:  %.4f, radius difference: %.4f \n'],...
        metrics_mdeit.centroid(1),metrics_mdeit.centroid(2),metrics_mdeit.radius,centroid_diff_mdeit,radius_diff_mdeit);
end




%% Functions
function show_fem_transparent_edges(img)

hh = show_fem(img);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

end

function metrics = compute_greit_metrics_2d(img,n_points)
    img.calc_colours.npoints = n_points;

    X_hat = calc_slices(img,[inf,inf,0]);

    Xq_hat = (abs(X_hat) >= 0.25*max(abs(X_hat(:))));
    
    % Keep only the largest connected blob. Avoid small artifacts
    CC = bwconncomp(Xq_hat);
    numPixels = cellfun(@numel, CC.PixelIdxList);

    [~, idx] = max(numPixels);     % largest component
    Xq_hat_clean = false(size(Xq_hat));
    Xq_hat_clean(CC.PixelIdxList{idx}) = true;

    stats = regionprops(Xq_hat_clean, 'Centroid');

    centroid_pixels = stats.Centroid;
    
    figure
    hold on
    imagesc(Xq_hat_clean);
    plot(centroid_pixels(1),centroid_pixels(2),'rx');
    xlim([1,n_points]);ylim([1,n_points])
    hold off;

   function physical_length = pixel_length_to_physical(img,pixel_length)

        physical_length = max(img.fwd_model.nodes(:,1))-min(img.fwd_model.nodes(:,1));
        physical_height = max(img.fwd_model.nodes(:,2))-min(img.fwd_model.nodes(:,2));

        assert(abs(physical_length - physical_height)<1e-12); % assume square

        physical_length = pixel_length(1)*physical_length/n_points;
    end

    function physical_position = pixel_position_to_physical(img,pixel_position)
        
        assert(all(size(pixel_position) == [1,2]) | all(size(pixel_position) == [2,1]));

        physical_min_x = min(img.fwd_model.nodes(:,1));
        physical_min_y = min(img.fwd_model.nodes(:,2));
        
        physical_length = max(img.fwd_model.nodes(:,1))-min(img.fwd_model.nodes(:,1));
        physical_height = max(img.fwd_model.nodes(:,2))-min(img.fwd_model.nodes(:,2));
        
        assert(abs(physical_length - physical_height)<1e-12); % assume square

        pos_x = physical_min_x + pixel_position(1)*physical_length/n_points;
        pos_y = physical_min_y + pixel_position(2)*physical_height/n_points;
        
        physical_position = [pos_x,pos_y];
    end
    
    % Physical centroid
    centroid = pixel_position_to_physical(img,centroid_pixels);
    metrics.centroid = centroid;

    %% Radius = Mean distance from centroid to boundary pixels
     
    % First need to find the boundary pixels
    [boundary_pixels_y, boundary_pixels_x] = find(bwperim(Xq_hat_clean)); % [y,x] is correct, not [x,y]
    hold on
    plot(boundary_pixels_x, boundary_pixels_y, 'r.', 'MarkerSize', 8)
    hold off

    radius = mean(sqrt(...
        (boundary_pixels_x - centroid_pixels(1)).^2 + ...
        (boundary_pixels_y - centroid_pixels(2)).^2));
    
    temp = pixel_length_to_physical(img,[radius,0]);
    physical_radius = temp;

    metrics.radius = physical_radius;
end



function ensure_spm12()
    % DESCRIPTION: 
    % Ensure SPM12 is available locally in the external folder.Prompt the
    % user to clone it if not.

    repoURL = 'https://github.com/spm/spm12.git';
    localDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                        '..', 'external', 'spm12');

    if isfolder(localDir)
        fprintf('SPM12 found at:\n  %s\n\n', localDir);

        fprintf('Added SPM12 to path.\n');
        addpath(localDir);

        check_dependencies_spm();

        return
    end

    fprintf('SPM12 not found at:\n  %s\n\n', localDir);
    reply = input('Do you want to clone the SPM12 repository now? [y/N]: ','s');

    if ~strcmpi(reply,'y')
        error('SPM12 is required to run this script.');
    end

    % Check that git exists
    [status,~] = system('git --version');
    if status ~= 0
        error('Git is not available on this system.');
    end

    % Create external folder if missing
    parentDir = fileparts(localDir);
    if ~isfolder(parentDir)
        mkdir(parentDir);
    end

    fprintf('Cloning SPM12...\n');
    cmd = sprintf('git clone %s "%s"', repoURL, localDir);
    status = system(cmd);

    if status ~= 0
        error('Failed to clone SPM12.');
    end

    fprintf('SPM12 successfully cloned.\n');

    fprintf('Add SPM12 to path.\n');
    addpath(localDir);

end


function ensure_noise_sample_file()
    % DESCRIPTION: 
    % Check if EnvironmentalNoise folder exists (does not check if it has the
    % necessary data). If it doesn't, propmt the user to download it
    
    script_folder = fileparts(mfilename('fullpath'));
        
    noise_dir = fullfile(script_folder,'EnvironmentalNoise');

    if isfolder(noise_dir)
        fprintf('EnvironmentalNoise folder found at :\n  %s\n\n', noise_dir);
        addpath(genpath(noise_dir))
        return;
    end
    
    noise_zip_url = 'https://zenodo.org/records/7872660/files/EnvironmentalNoise.zip';
    zip_file = fullfile(script_folder, 'EnvironmentalNoise.zip');
    
    fprintf('EnvironmentalNoise data not found.\n');
    reply = input('Do you want to download it now (~4Gb)? [y/N]: ', 's');
    
    if ~strcmpi(reply,'y')
        error('EnvironmentalNoise data is required to run this script.');
    end

    fprintf('Downloading EnvironmentalNoise.zip...\n');
    try
        websave(zip_file, noise_zip_url);
    catch ME
        error('Failed to download EnvironmentalNoise.zip:\n%s', ME.message);
    end

    fprintf('Unzipping...\n');
    try
        unzip(zip_file, script_folder);
    catch ME
        delete(zip_file);
        error('Failed to unzip EnvironmentalNoise.zip:\n%s', ME.message);
    end

    delete(zip_file);

    if ~isfolder(noise_dir)
        error('EnvironmentalNoise folder not found after unzip.');
    end

    fprintf('EnvironmentalNoise successfully installed at:\n  %s\n\n', noise_dir);
    addpath(genpath(noise_dir))
end


function check_dependencies_spm()
    
    return;
end