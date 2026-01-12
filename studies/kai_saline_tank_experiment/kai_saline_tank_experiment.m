clc; clear all; close all;

%% Prepare workspace
% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

% Have to add the functions path manually so prepare_workspace runs
grandparent_folder =fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

data_folder = strcat(script_folder ,'\data');

clc;
rng(1)

% PROBLEM:

% GCV results are different depending on the mesh used for the forward
% model resulting from mesh convergence. Verify if the forward magnetic
% field computed with this mesh is actually converging. Might be a mistake
% in the R matrices maybe? 

% Additionally. If I erase the models from the model folder, visually the
% results are different and also the optimal lambda from GCV. When I erase,
% the optimal lambda is 1.1e-8, and when I run without erasing its 4e-10. 
% IT SURELLY HAS SOMETHING TO DO WITH LOADING THE MODEL FROM MEMORY!!!!!!

% The problem might be the backslash inversion in the GCV code. Use pcg
% instead? Doubt this. Its probably what is in CAPS LOCK above.



% CODE in need of testing:

% -> Cylinder perturbation sometimes gets NETGEN stuck. This is a problem
% in EIDORS ng_mk_cyl_models when using the "extra" input argument for a
% cylinder and orthobrick.

% TODO/CORRECTIONS:

% -> Add a failure state to the inverse solver, that fails when pcg fails. It
% seems that sometimes the inverse solve goes by very quick for the pcg to
% be working correctly

% -> Check if pcg number of iterations/tolerance is too big to solve
% correctly for noise

% -> Re-scale conductivity so conductivity of ping pong ball is bigger than
% machine precision

% -> When no noise is given, EIDORS solver for EIT complains that the GCV
% system is ill-conditioned. The results come out meaningless for both EIT
% and MDEIT in that case

% -> If user gives sensor locations, provide check that these do not
% overlap

% -> Should have fwd_model parameters and img parameters separatly,
% todo.

% -> Run a mesh convergence study. Almost DONE, but need to implement
% convergence criterion.

% ->  For now, can't seem to make EIDORS want to place circular electrodes,
% so they are square currently.

% -> Gotta figure out, when building a mesh with a material, how to assign
% it a name

% -> if the Jacobian has been computed for some reason, store it in the
% forward model instead of recomputing, but check if size is reasonable
% before.

% -> Look into tensorprod function in matlab to accelerate Jacobian row
% assembly

% -> It seems that as I add more elements, the solution for the 1-axis
% MDEIT in 2D becomes messy. I think there might be something wrong with
% the 2D solver perhaps.

% -> electrodeRadius property is actually the side length ( or diameter),
% correct that


% NOTES:

% The electrodes used in the tank experiment are only mentioned to be
% Ag/AgCl EEG electrodes. I suppose that makes the circular electrodes, but
% the radius is unknow, and depends on the manufacterer apparently, but
% from what I've seen on the internet, it's around the range of {4,8,12}mm.
% Lets consider 8mm.

% Requesting circular electrodes makes EIDORS break for some reason. Look
% into that later.

% 0 S/m conductivity perturbation makes EIDORS forward solution break. Lets
% use the values from this website: https://matmake.com/properties/electrical-conductivity-of-polymers-and-plastics.html
% for Acrylonitrile butadiene styrene (ping pong ball plastic) 10e-13 S/m =
% 1e-15 S/mm = 1e-12 mS/mm

% Can't easily find the electrical conductivity of the aqueous NaCl
% solution used, so I checked this site: https://diverdi.colostate.edu/all_courses/CRC%20reference%20data/electrical%20conductivity%20of%20aqueous%20solutions.pdf
% and take the conductivity of NaCl (aq) at 0.2% mass to be 2/5 of the
% conductivity for the 0.5% mass case: 8.2*2/5 = 3.28 mS/cm = 3.28e-4 S/mm
% = 3.28e-1 mS/mm

% When the plastic conductivity is very small, the forward problem matrix
% is ill-conditioned, even for EIDORS, so I changed to values of
% conductivity to mS/mm units (THIS IS NOT TRUE, LOOK at next note). In these units, current density is given in
% A/m^2 if the electric field is given in V/m.

% Turns out that the problem was that EIDORS uses the hp^2*RtR in the
% normal equations, while the generalized_cross_validation function
% consideres hp*RtR. So the hyperparameter that I was feeding to EIDORS was
% very small which causes J'*J+hp^2*RtR to not be positive definite and the
% EIDORS solver using linsolve would error. I kept the changes in units
% anyway. The correction is giving the EIDORS inverse solver sqrt(hp)

% I was using the experiment parameters from the two dimensional tank
% experiment, but later in his thesis, Kai has a two dimensional simulated
% experiment. This is better to compare his results to mine. I had to
% change the data from the experiment to the one found on page 248. Some of
% the data for this experiment is missing tho, so this script is a mix
% between both his real and simulated tank experiment.

% I changed how I add noise, to add directly to the difference signal
% instead to the absolute measured magnetic fields and voltages. This is
% not very realistic, but is better to compare EIT and MDEIT on an equal
% footing. 



%Problems: Number of magnetic sensors in 2D computational experiment is not
%mentioned

%The electrode contact impedance used in the FEM model is not mentioned.
%In CEM article, they use 58.0,35.0,15.0,7.5 Ohm cm^2 (before discussion section),
%which would be 58*(0.01)^2 = 0.0058 Ohm m^2. With these value, the contact
%impedance for 8mm radius electrodes is around 0.0058/(2*pi*8e-3^2) = 14.5 Ohm. In EIDORS, they
%mention that the contact impedance is 20 Ohm, so it seems consistent with
%that

% The noise I was adding to EIT was not that comparable to the noise in
% MDEIT, because the signal is the difference between inhomogeneous and
% homogeneous measurements, and that difference becomes comparable to the
% noise quicker for MDEIT than for EIT while reducing the SNRdb. For now,
% unrealistically, add noise to the difference signal directly.

%For some reason, Kai's noise correction works if I divide by the std the
%elem_data-mean(elem_data) instead of the elem_data directly. AH!!!! It's
%because the GN method is outputting not the difference data, but the
%homogeneous solution plus the difference!!!!!!!!!!!!!!!!!!!! GEEEZE. Its
%correct now

%TODO: The reconstruction changes a bit with increasing mesh refinement, or maybe
%its just the effect of computing the R matrices with quadrature and not
%with Matlab's integral2. Check if reconstructions with coarser meshes are
%the same with the quadrature R matrices!

% Question: Not sure if it makes sense to use a difference reconstruction
% method when the background conductivity and anomaly conductivity are so
% different. We are not in the regime where \Delta\sigma is small. Why did
% Kai not consider absolute reconstruction?


%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
%% Assign the parameters for several models (should create utility functions for this)

% For small SNR, noise level might be bigger in magnitude than difference
% data. Does that measurement even make sense?
SNRdb = 20;


% minsz_mesh_convergence = 0.05;
% maxsz_mesh_convergence = 0.2;
% num_meshes_mesh_convergence = 5;

% maxsz_reconstruction = 0.05;

minsz_mesh_convergence = 0.05;
maxsz_mesh_convergence = 0.05;
num_meshes_mesh_convergence = 1;

maxsz_reconstruction = 0.1;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

%% Define forward model (2D real tank experiment)

% maxsz is determined by mesh convergence study
model_parameters.maxsz = 10e-3/l0;

% This data came from Kai's Thesis
model_parameters.is2D = true;

% Tank dimensions in characteristic units
model_parameters.isCylindrical = true; % page 163
model_parameters.height=70e-3/l0; %(m) page 163 
model_parameters.radius=40e-3/l0; %(m) page 163 

%Electrode configuration and properties
model_parameters.numOfRings = 1;
model_parameters.numOfElectrodesPerRing = 16; 

% what's the size of the electrodes? He mentions EEG for the tank experiment, but what for the simulated experiment?
model_parameters.electrodeRadius = 8e-3/l0; %(mm) unknown from the thesis, but good estimate for Ag/AgCl EEG electrode 
assert(model_parameters.electrodeRadius*model_parameters.numOfElectrodesPerRing < 2*pi*model_parameters.radius)
model_parameters.electrodeContactImpedance = 0.0058/z0; 

% Magnetic sensor configuration and properties
model_parameters.numOfSensors = 25; %(mm) pg 172
model_parameters.sensorRadius = 70e-3/l0;  %(mm) pg 172
model_parameters.mu0 = 1; %With this definition, we're computing the magnetic field intensity H in units of I0/l0;


anomaly_conductivity = 1e-12/sigma0; % (0 according to Kai's thesis, but check notes)
anomaly_position = [20,0,35]*1e-3/l0; %  [0,0,35]*1e-3/l0; % pg 172
anomaly_radius = 25/2*1e-3/l0; % 25mm of diameter, pg 172

% This is deprecated, the .material property overrides this behaviour,
% except for the anomaly conductivity
model_parameters.anomaly = struct(...
    'type','spherical',...
    'conductivity',anomaly_conductivity,...
    'radius',anomaly_radius,...
    'position',anomaly_position);

% A material defines the geometry of that material inside the domain, not
% its properties, like conductivity.
model_parameters.material = struct(...
    'type','cylindrical',...
    'name','plastic_cylinder',...
    'radius',anomaly_radius,...
    'position',anomaly_position);

% electrode height 35 mm is correct given that the ring is placed in the
% middle of the height of the cylinder even though there is no parameter of
% model_parameters controlling it

% Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
options = {};
%% Mesh convergence study to obtain the "correct" forward model

[fmdl_converged,u_norms,b_norms,n_elems,n_nodes] = ...
    mesh_convergence(model_parameters,model_folder,minsz_mesh_convergence,maxsz_mesh_convergence,num_meshes_mesh_convergence,'multiple_stimulation');

fmdl = fmdl_converged;

dim = size(fmdl.nodes,2);

%% Assign the stimulation pattern to this forward model

stimulation = mk_stim_patterns(numel(fmdl.electrode),1,inj,meas,options,current_amplitude);
fmdl.stimulation = stimulation;
%% Make images from forward models and forward solve

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);

% Add plastic cylinder
imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

%% Generate data from these forward models, but consider a different one for reconstruction

% Forward solve
[datah,uh] = fwd_solve_mdeit(imgh);
Bzh = datah.Bz(:);

% Forward solve
[datai,ui] = fwd_solve_mdeit(imgi);
Bzi = datai.Bz(:);

Bzh_real = Bzh;
Bzi_real = Bzi;

ui_real = ui;
uh_real = uh;

%% Visualize pattern 1 and conductivity anomaly

% Conductivity

% Electric potential
img_visualize_u = rmfield(imgi, 'elem_data');
img_visualize_u.node_data = ui.volt(:,1);

img_visualize_J = rmfield(imgi, 'elem_data');

% Current density norm
Jx = -imgi.elem_data.*imgi.fwd_model.G.Gx*ui.volt(:,1);
Jy = -imgi.elem_data.*imgi.fwd_model.G.Gy*ui.volt(:,1);
Jz = -imgi.elem_data.*imgi.fwd_model.G.Gz*ui.volt(:,1);

normJ = sqrt(Jx.^2+Jy.^2+Jz.^2);

img_visualize_J.elem_data = normJ;

% Set up colorbars

% Compute colors and scaling limits
imgi.calc_colours.ref_level = background_conductivity;  %  centre of the colour scale
imgi.calc_colours.clim = abs(background_conductivity-anomaly_conductivity); %  max diff from ref_level

img_visualize_u.calc_colours.ref_level = mean(ui.volt(:,1));
img_visualize_u.calc_colours.clim = max(ui.volt(:,1)-mean(ui.volt(:,1)));

img_visualize_J.calc_colours.ref_level = mean(normJ);
img_visualize_J.calc_colours.clim = max(normJ);

% calc_colours
figure('Position',[200,100,1200,600]);
hold on
subplot(2,2,1)

hh = show_fem(imgi);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

plot_sensors(fmdl_converged)
if dim ==2
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius])
    axis square
else
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius...
        0 model_parameters.height])
end
eidors_colourbar(imgi);
title('$\sigma_{ground \: truth}$ - units $z_0^{-1}l_0$','Interpreter','latex')
box on;
subplot(2,2,2)
hold on

hh = show_fem(img_visualize_u);              
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

plot_sensors(fmdl_converged)
if dim ==2
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius])
    axis square
else
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius...
        0 model_parameters.height])
end
hold off
eidors_colourbar(img_visualize_u);
box on;
title('$u$ - injection pattern $1$ - units $z_0 l_0^{-2} I_0$','Interpreter','latex')
hold off

subplot(2,2,3)
hold on

hh = show_fem(img_visualize_J);            
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

plot_sensors(fmdl_converged)
if dim ==2
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius])
    axis square
else
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius...
        0 model_parameters.height])
end
cb = eidors_colourbar(img_visualize_J);

hold off
box on;
title('$||\vec{J}||$ - injection pattern $1$ - units $I_0l_0^{-2}$','Interpreter','latex')
hold off

subplot(2,2,4)
hold on
plot(1:numel(Bzh),Bzh);

xlabel('measurement index','Interpreter','latex');
ylabel('$H^z$','Interpreter','latex');

hold off
box on;
title('$H^z_{homogeneous}$ - units $I_0l_0^{-1}$','Interpreter','latex')
hold off
axis square

pause(1e-5);

%% Add measurement noise
% Bzh = add_measurement_noise(Bzh,SNRdb);
% Bzi = add_measurement_noise(Bzi,SNRdb);

% uh.meas = add_measurement_noise(uh.meas,SNRdb);
% ui.meas= add_measurement_noise(ui.meas,SNRdb);

dB = add_measurement_noise_difference(Bzi,Bzh,SNRdb);
du = add_measurement_noise_difference(ui.meas,uh.meas,SNRdb);

Bzi = Bzh+dB;
ui.meas = uh.meas+du;

%% Plot real data vs noisy data

dB = Bzi-Bzh;
du = ui.meas-uh.meas;

dB_real = Bzi_real-Bzh_real;
du_real = ui_real.meas-uh_real.meas;

noise_mdeit = dB-dB_real;
noise_eit = du-du_real;

figure
subplot(3,2,1)
hold on
plot(Bzi_real,'b');
plot(Bzh_real,'g-.');
plot(Bzi,'r.');

ylabel('$H^z_i$','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
legend('Real inh','Real h','Noisy inh')
title('$H^z_i$','Interpreter','latex');

subplot(3,2,2)
hold on
plot(ui_real.meas,'b');
plot(uh_real.meas,'g-.');

plot(ui.meas,'r.');
ylabel('$u_i$','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
legend('Real inh','Real h','Noisy inh')
title('$u_i$','Interpreter','latex');


subplot(3,2,3)
hold on
plot(dB_real,'b');
plot(dB,'r.');

id_mdeit = find(abs(noise_mdeit) == max(abs(noise_mdeit)));
line = [id_mdeit,dB_real(id_mdeit);id_mdeit,dB(id_mdeit)];
plot(line(:,1),line(:,2),'k-','LineWidth',1,'LineStyle','--');
plot(line(1,1),line(1,2),'k*');

ylabel('$\Delta H^z$','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
legend('Real','Noisy','Max Noise')
title('$H^z_i-H^z_h$','Interpreter','latex');

subplot(3,2,4)
hold on
plot(du_real,'b');
plot(du,'r.');

id_eit = find(abs(noise_eit) == max(abs(noise_eit)));
line = [id_eit,du_real(id_eit);id_eit,du(id_eit)];
plot(line(:,1),line(:,2),'k-','LineWidth',1,'LineStyle','--');
plot(line(1,1),line(1,2),'k*');


ylabel('$\Delta u$','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
legend('Real','Noisy','Max Noise')

title('$u_i-u_h$','Interpreter','latex');

subplot(3,2,5)
hold on
plot(noise_mdeit,'b');

ylabel('noise','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
title('noise','Interpreter','latex');

subplot(3,2,6)
hold on
plot(noise_eit,'b');

ylabel('noise','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
title('noise','Interpreter','latex');

relative_change_eit = abs(du_real)./(abs(uh_real.meas));
relative_change_mdeit = abs(dB_real)./(abs(Bzh_real));

% Print some diagonostics
fprintf('------------------------------------------\n')
max_noise_eit = abs(noise_eit(id_eit));
max_noise_mdeit = abs(noise_mdeit(id_mdeit));

fprintf('Max noise EIT [%%]: %.3g %% of max signal\n',...
    100*max_noise_eit/max(du_real));
fprintf('Max noise MDEIT [%%]: %.3g %% of max signal\n',...
    100*max_noise_mdeit/max(dB_real));


pause(1e-5)

%% Emit a warning when the noise is of comparable size as the difference signal

% amplitude_mdeit = max(abs(Bzi-mean(Bzi)));
% signal = Bzi_real-Bzh_real;
% signal_amplitude_mdeit = max(abs(signal-mean(signal)));
% noise_amplitude_mdeit = 2*amplitude_mdeit/10^(SNRdb/20); % times 2 because in the worst case the noise from Bzh and Bzi adds up
% 
% if noise_amplitude_mdeit>signal_amplitude_mdeit
% 
%     warning('MDEIT noise amplitude is bigger than signal')
%     choice = questdlg('MDEIT noise amplitude is bigger than signal. Proceed with computation?', ...
%         'Confirmation', ...
%         'Yes','No','No');
%     switch choice
%         case 'Yes'
%             %continue
%         case 'No'
%             error('User aborted.');
%     end
% 
% end
% 
% amplitude_eit = max(abs(ui.meas-mean(ui.meas)));
% signal = ui_real.meas-uh_real.meas;
% signal_amplitude_eit = max(abs(signal-mean(signal)));
% noise_amplitude_eit = 2*amplitude_eit/10^(SNRdb/20); % times 2 because in the worst case the noise from Bzh and Bzi adds up
% 
% if noise_amplitude_eit>signal_amplitude_eit
% warning('EIT noise amplitude is bigger than signal')
%     choice = questdlg('EIT noise amplitude is bigger than signal. Proceed with computation?', ...
%         'Confirmation', ...
%         'Yes','No','No');
%     switch choice
%         case 'Yes'
%             %continue
%         case 'No'
%             error('User aborted.');
%     end
% end

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

% 3-axis imdl
imdl_mdeit_3 = rmfield(imdl_mdeit_1, 'select_sensor_axis');
imdl_mdeit_3.recon_mode = 'mdeit3';

% eit imdl
% eidors_cache('clear_all'); %For safety, because EIDORS is using the cache forward solve from a different forward model.
imdl_eit = mk_common_model('a2c2',8); % Will replace most fields
imdl_eit.jacobian_bkgnd = struct('value',background_conductivity); %If I
% use this, the solution blows up
imdl_eit.fwd_model = fmdl_reconstruction; %Use a different forward model for reconstruction
imdl_eit.recon_mode = 'eit';
imdl_eit.reconst_type = 'difference';
imdl_eit.RtR_prior = @prior_tikhonov;

%% Find optimal regularization parameter with GCV
lambda_vector = logspace(-15,3,30);

[lambda_optimal_mdeit1,optimal_id_mdeit_1,V_mu_mdeit1,~] = generalized_cross_validation(imdl_mdeit_1,Bzi-Bzh,lambda_vector);

imdl_mdeit_1.hyperparameter = struct('value',lambda_optimal_mdeit1);

n_eit = numel(fmdl.stimulation)*size(fmdl.stimulation(1).meas_pattern);

[lambda_optimal_eit,optimal_id_eit,V_mu_eit,~] = generalized_cross_validation(imdl_eit,ui.meas-uh.meas,lambda_vector);

% Since EIDORS uses hp^2*RtR, but the generalized_cross_validation function
% uses hp*RtR, we have to correct for that.
imdl_eit.hyperparameter = struct('value',sqrt(lambda_optimal_eit));

%% Show the generalized cross validation results
figure;
subplot(1,2,1);
hold on
plot(lambda_vector,V_mu_eit);
plot(lambda_vector(optimal_id_eit),V_mu_eit(optimal_id_eit),'r.');
box on;grid on;grid minor;
set(gca,'YScale','log');
set(gca,'XScale','log');
xlim([min(lambda_vector),max(lambda_vector)])

xlabel('$\lambda$','Interpreter','latex')
ylabel('$V(\lambda)$','Interpreter','latex')

title('GCV EIT','Interpreter','latex')

subplot(1,2,2);
hold on
plot(lambda_vector,V_mu_mdeit1);
plot(lambda_vector(optimal_id_mdeit_1),V_mu_mdeit1(optimal_id_mdeit_1),'r.');
set(gca,'YScale','log');
set(gca,'XScale','log');
box on;grid on;grid minor;
xlim([min(lambda_vector),max(lambda_vector)])
xlabel('$\lambda$','Interpreter','latex')
ylabel('$V(\lambda)$','Interpreter','latex')
title('GCV MDEIT','Interpreter','latex')

%% Perform reconstruction

% Reconstruct (difference) for 1-MDEIT
img_output_mdeit_1 = inv_solve_mdeit(imdl_mdeit_1,Bzh,Bzi);

% Reconstruct (difference) for EIT
img_output_eit = inv_solve(imdl_eit,uh,ui);

%% Apply KAI's noise correction

close all
num_noise_repetitions = 100;

% Solve Tykhonov regularized normal equations using SVD for noisy inputs

if isfield(img_output_mdeit_1,'jacobian')
    J_mdeit = img_output_mdeit_1.jacobian;
else
    error('Expected to have jacobian here');
end

img_eit = mk_image(fmdl_reconstruction,background_conductivity);

J_eit = calc_jacobian(img_eit);

[U_eit,S_eit,V_eit] = svd(J_eit,'econ');
[U_mdeit,S_mdeit,V_mdeit] = svd(J_mdeit,'econ');

M_mdeit = zeros(n_elem,num_noise_repetitions);
M_eit = zeros(n_elem,num_noise_repetitions);

figure; close all;

noise_level_mdeit= max(abs(dB-mean(dB)))/10^(SNRdb/20);
noise_level_eit = max(abs(du-mean(du)))/10^(SNRdb/20);

for t = 1:num_noise_repetitions
    
    fprintf('Solving for noise realization %i\n',t);

    data_noisy_mdeit = noise_level_mdeit*randn(size(Bzi,1),1);
    data_noisy_eit = noise_level_eit*randn(size(ui.meas,1),1);
    
    % z_mdeit = (S_mdeit'*S_mdeit+imdl_mdeit_1.hyperparameter.value*speye(size(S_mdeit'*S_mdeit)))\...
    %     (S_mdeit'*U_mdeit'*data_noisy_mdeit);
    % sigma_mdeit_noisy = V_mdeit*z_mdeit;
    % 
    % z_eit = (S_eit'*S_eit+imdl_eit.hyperparameter.value^2*speye(size(S_eit'*S_eit)))\...
    %     (S_eit'*U_eit'*data_noisy_eit);
    % sigma_eit_noisy = V_eit*z_eit;
    
    sv_mdeit = diag(S_mdeit)+imdl_mdeit_1.hyperparameter.value./diag(S_mdeit);
    sigma_mdeit_noisy = V_mdeit*diag(1./sv_mdeit)*U_mdeit'*data_noisy_mdeit;
    
    % sigma_mdeit_noisy = ...
    %     V_mdeit*...
    %     inv(S_mdeit*S_mdeit+imdl_mdeit_1.hyperparameter.value*eye(size(S_mdeit)))*...
    %     S_mdeit*U_mdeit'*data_noisy_mdeit; %this approach is the same as above
    
    % sigma_mdeit_noisy_2 = linsolve(...
    %     J_mdeit'*J_mdeit+imdl_mdeit_1.hyperparameter.value*speye(n_elem),...
    %     J_mdeit'*data_noisy_mdeit);%this approach is the same as above
    
    % Careful here, have to consider that for EIT, the solver solved
    % J'J+hp^2*RtR and not J'J+hp*RtR
    sv_eit = diag(S_eit)+imdl_eit.hyperparameter.value.^2./diag(S_eit);
    sigma_eit_noisy = V_eit*diag(1./sv_eit)*U_eit'*data_noisy_eit;

    % sigma_eit_noisy_2 = inv_solve(imdl_eit,zeros(size(data_noisy_eit)),data_noisy_eit);
    % sigma_eit_noisy_2 =  sigma_eit_noisy_2.elem_data; %this approach is the same as above

    M_mdeit(:,t) = sigma_mdeit_noisy;
    M_eit(:,t) = sigma_eit_noisy;
    
    % Debug
    % sigma_std_deviation_mdeit = std(M_mdeit(:,1:t),[],2);
    % sigma_std_deviation_eit = std(M_eit(:,1:t),[],2);
    % 
    % sigma_std_deviation_eit = sigma_std_deviation_eit/max(sigma_std_deviation_eit);
    % sigma_std_deviation_mdeit  = sigma_std_deviation_mdeit/max(sigma_std_deviation_mdeit);
    % 
    % 
    % subplot(1,2,1);
    % cla;
    % plot(sigma_std_deviation_eit);
    % subplot(1,2,2);
    % cla;
    % plot(sigma_std_deviation_mdeit);
    % 
    % pause(1e-5)

end

sigma_std_deviation_mdeit = std(M_mdeit,[],2);
sigma_std_deviation_eit = std(M_eit,[],2);

img_output_eit_corrected = img_output_eit;
img_output_eit_corrected.elem_data = img_output_eit.elem_data./sigma_std_deviation_eit;

img_output_mdeit_1_corrected = img_output_mdeit_1;
img_output_mdeit_1_corrected.elem_data = img_output_mdeit_1.elem_data./sigma_std_deviation_mdeit;


%% Try Kaipio's approach of computing posterior covariance matrix

Gamma_pr = background_conductivity*diag(ones(n_elem,1));
Gamma_noise = noise_level_mdeit*eye(size(J_mdeit,1));

Gamma_post = inv(J_mdeit'*inv(Gamma_noise)*J_mdeit+inv(Gamma_pr)); %Moore Penrose pseudo inverse
% Gamma_post should be symmetric, so force symmetry. This stops problems
% with eig

Gamma_post = (Gamma_post+Gamma_post')/2;

img_output_mdeit_1_corrected_3 = img_output_mdeit_1;
img_output_mdeit_1_corrected_3.elem_data = img_output_mdeit_1_corrected_3.elem_data./sqrt(diag(Gamma_post));

% img_output_mdeit_1_corrected_3 = img_output_mdeit_1;
% img_output_mdeit_1_corrected_3.elem_data = img_output_mdeit_1_corrected_3.elem_data./sqrt(diag(Gamma_post));

img_output_mdeit_1_corrected_4 = img_output_mdeit_1;

[V1,D1] = eig(Gamma_post);
gamma = 0.9*max(diag(D1));
g = diag(D1)./(diag(D1) + gamma);

img_output_mdeit_1_corrected_4.elem_data =...
    V1 * diag(g) * (V1' * img_output_mdeit_1_corrected_4.elem_data) + ...   %projection into eigenbasis (shrink high-variance modes by factor of g)
    (eye(size(V1*V1')) - V1*V1')*img_output_mdeit_1_corrected_4.elem_data;  %orthogonal projection into kernel (don't do anything)

%% Save data

file_name = strcat(script_folder,'/data/kai_saline_tank_SNR_',num2str(SNRdb));

save(file_name,"imgi",...
    "imdl_eit","imdl_mdeit_1",...
    "img_output_mdeit_1","img_output_eit",...
    "sigma_std_deviation_eit","sigma_std_deviation_mdeit");

%% Show reconstruction
figure('Position',[200,200,1000,500]);

subplot(2,2,1)
show_fem(imgi)
box on
title('Ground-Truth','Interpreter','latex')
subplot(2,2,2)
show_fem(img_output_eit)
eidors_colourbar(img_output_eit)
box on
title('EIT','Interpreter','latex')
subplot(2,2,3)
show_fem(img_output_mdeit_1)
eidors_colourbar(img_output_mdeit_1)
box on
title('$1$-MDEIT','Interpreter','latex')
subplot(2,2,4)
if dim == 3
    show_fem(img_output_mdeit_3)
end
box on
title('$3$-MDEIT','Interpreter','latex')

%% Show reconstruction (noise corrected)

figure('Position',[100,200,1300,500]);

subplot(2,3,1)
show_fem(img_output_eit_corrected)
eidors_colourbar(img_output_eit_corrected)
box on
title('EIT (Noise corrected)','Interpreter','latex')

subplot(2,3,2)
show_fem(img_output_mdeit_1_corrected)
box on
title('$1$-MDEIT (Noise corrected - Kai TSVD)','Interpreter','latex')

subplot(2,3,3)
imgr = rmfield(img_output_eit,'elem_data');
data = img_output_mdeit_1.elem_data.^(1);

% data = data./sigma_std_deviation_eit;
imgr.elem_data = data;
show_fem(imgr)
eidors_colourbar(imgr)
box on
title('$1$-MDEIT (Noise corrected using EIT std)','Interpreter','latex')

% subplot(2,3,4)
% show_fem(img_output_mdeit_1_corrected_3)
% box on
% title('$1$-MDEIT (Noise corrected - Kaipio, std from covariance matrix)','Interpreter','latex')

subplot(2,3,4)
imgr = rmfield(img_output_eit,'elem_data');
imgr.elem_data = img_output_eit.elem_data./sigma_std_deviation_mdeit;
show_fem(imgr)
box on
title('EIT (Noise corrected using MDEIT std)','Interpreter','latex')

% subplot(2,3,5)
% show_fem(img_output_mdeit_1_corrected_4)
% box on
% title('$1$-MDEIT (Noise corrected - Kaipio, kill high variance modes)','Interpreter','latex')

subplot(2,3,5)
imgr = rmfield(img_output_eit,'elem_data');
imgr.elem_data = sigma_std_deviation_eit;
show_fem(imgr)
eidors_colourbar(imgr)
box on
title('EIT (Kai std)','Interpreter','latex')

subplot(2,3,6)
imgr = rmfield(img_output_mdeit_1,'elem_data');
imgr.elem_data = sigma_std_deviation_mdeit;
show_fem(imgr)
eidors_colourbar(imgr)
box on
title('$1$-MDEIT (Kai std)','Interpreter','latex')
%% FUNCTIONS

function data_noisy = add_measurement_noise(data,SNRdb)
assert(isnumeric(data) && isvector(data),'data must be a numerical vector');
assert(size(data,1)>=size(data,2),'data must be a column vector');

if isempty(SNRdb)
    data_noisy = data;
    return;
else
    data_amplitude = std(data);
    noise_level= data_amplitude/10^(SNRdb/20);
    data_noisy = data+noise_level*randn(size(data));
end
end


function data_noisy = add_measurement_noise_difference(datai,datah,SNRdb)

if isempty(SNRdb)
    data_noisy = datai-datah;
    return;
else
    d = datai-datah;
    d_amplitude = max(abs(d-mean(d)));

    noise_level= d_amplitude/10^(SNRdb/20);

    data_noisy = d+noise_level*randn(size(d));
end
end



%%
figure
hh = show_fem(img_output_mdeit_1_corrected);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);
axis off