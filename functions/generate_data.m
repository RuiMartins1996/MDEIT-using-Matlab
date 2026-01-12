function  [out_data_b,out_data_u,snr,noise_b,noise_u] = generate_data(options,img1,varargin)

assert(strcmp(img1.type,'image'),'First argument must be image');
assert(nargin<4,'Wrong number of arguments');

assert(isfield(options,'noise_structure'));

if nargin==2
    recon_mode = 'absolute';
elseif nargin == 3
    recon_mode = 'difference';
    img2 =  varargin{1};
end

[noise_structure] = validate_options(options);

% Get number of measurements per stimulation
num_of_measurements = zeros(1,numel(img1.fwd_model.stimulation));
for i = 1:numel(img1.fwd_model.stimulation)
    num_of_measurements(i) = size(img1.fwd_model.stimulation(i).meas_pattern,1);
end

num_of_injections = numel(img1.fwd_model.stimulation);

% Number of measurements per stimulation must be the same for
% this to work
if not(all(num_of_measurements == num_of_measurements(1)))
    error('Number of measurements per stimulation must be the same');
end

num_of_eit_channels = num_of_measurements(1);
num_of_mdeit_channels_per_dim = numel(img1.fwd_model.sensors);

if isfield(noise_structure,'snr')
    snr = noise_structure.snr;
else
    snr = 0;
end

switch recon_mode
    case 'difference'
        % Get data
        [datah_b,datah_u] = fwd_solve_mdeit(img1);
        [datai_b,datai_u] = fwd_solve_mdeit(img2);

        dBx = datai_b.Bx(:)-datah_b.Bx(:);
        dBy = datai_b.By(:)-datah_b.By(:);
        dBz = datai_b.Bz(:)-datah_b.Bz(:);

        du = datai_u.meas-datah_u.meas;

        if any(strcmp(noise_structure.type,{'white','pink','brown'}))

            assert(isfield(noise_structure,'sampling_rate'),'must have field sampling_rate here');
            n_samples_epoch = noise_structure.epoch_time*noise_structure.sampling_rate; % number of samples in the epoch
            noise_structure.n_samples_epoch = n_samples_epoch;
            
            dBx_noisy = dBx;
            dBy_noisy = dBy;
            dBz_noisy = dBz;

            noise_amplitude_x = signal_amplitude(dBx) / 10^(snr/20);
            noise_amplitude_y = signal_amplitude(dBy) / 10^(snr/20);
            noise_amplitude_z = signal_amplitude(dBz) / 10^(snr/20);

            noise_structure.num_of_channels = num_of_mdeit_channels_per_dim;
            
            noise_x = zeros(num_of_mdeit_channels_per_dim*num_of_injections,noise_structure.n_samples_epoch);
            noise_y = zeros(num_of_mdeit_channels_per_dim*num_of_injections,noise_structure.n_samples_epoch);
            noise_z = zeros(num_of_mdeit_channels_per_dim*num_of_injections,noise_structure.n_samples_epoch);

            for i = 1:num_of_injections
                ids  = (i-1)*num_of_mdeit_channels_per_dim+1:(i-1)*num_of_mdeit_channels_per_dim+num_of_mdeit_channels_per_dim;

                noise_x(ids,:) = generate_noise(noise_structure);

                noise_y(ids,:) = generate_noise(noise_structure);

                noise_z(ids,:) = generate_noise(noise_structure);
               
            end
            
            % Scale noise such that SNR is exactly the requested
            noise_b_x = noise_amplitude_x*mean(noise_x,2)./signal_amplitude(mean(noise_x,2));
            noise_b_y = noise_amplitude_y*mean(noise_y,2)./signal_amplitude(mean(noise_y,2));
            noise_b_z = noise_amplitude_z*mean(noise_z,2)./signal_amplitude(mean(noise_z,2));

            dBx_noisy = dBx + noise_b_x;
            dBy_noisy = dBy + noise_b_y;
            dBz_noisy = dBz + noise_b_z;
            
            % Sanity check
            if abs(compute_snr(dBz,dBz_noisy-dBz)-snr)>1e-3
                error('Somethings wrong');
            end
                
            noise_structure.num_of_channels = num_of_eit_channels;
            noise_amplitude_eit = signal_amplitude(du) / 10^(snr/20);
            
            noise_u = zeros(num_of_eit_channels*num_of_injections,noise_structure.n_samples_epoch);
            for i = 1:num_of_injections
                ids  = (i-1)*num_of_eit_channels+1:(i-1)*num_of_eit_channels+num_of_eit_channels;

                noise_u(ids,:) = generate_noise(noise_structure);   
            end
            
            % Scale noise such that SNR is exactly the requested
            noise_u = noise_amplitude_eit*mean(noise_u,2)./signal_amplitude(mean(noise_u,2));
            du_noisy = du+noise_u;

            % Sanity check
            if abs(compute_snr(du,du_noisy-du)-snr)>1e-3
                error('Somethings wrong');
            end
            

        elseif strcmp(noise_structure.type,'opm')

            % Generate noise time series for multiple channels
            noise_structure.num_of_channels = num_of_mdeit_channels_per_dim;
            num_of_injections = numel(img1.fwd_model.stimulation);
            noise_structure.num_of_measurements = num_of_injections;


            % Noise time series is num_of_channels x epoch length x num_of_injections
            % Choose a starting id of time series of raw noise
            time = noise_structure.D.time();
            dt = time(2) - time(1);          % sampling interval
            n_samples_epoch = round(noise_structure.epoch_time / dt); % number of samples in the epoch

            t_id = randi(nsamples(noise_structure.D)-6*num_of_injections*n_samples_epoch);
            [noise_time_series_x_h,t_id] = generate_opm_noise(noise_structure,t_id);
            [noise_time_series_y_h,t_id] =generate_opm_noise(noise_structure,t_id);
            [noise_time_series_z_h,t_id] = generate_opm_noise(noise_structure,t_id);

            [noise_time_series_x_i,t_id] = generate_opm_noise(noise_structure,t_id);
            [noise_time_series_y_i,t_id] = generate_opm_noise(noise_structure,t_id);
            [noise_time_series_z_i,~] = generate_opm_noise(noise_structure,t_id);

            % Create signal time series
            Bxh_noisy = zeros(size(noise_time_series_z_h));
            Byh_noisy = zeros(size(noise_time_series_z_h));
            Bzh_noisy = zeros(size(noise_time_series_z_h));

            Bxi_noisy = zeros(size(noise_time_series_z_h));
            Byi_noisy = zeros(size(noise_time_series_z_h));
            Bzi_noisy = zeros(size(noise_time_series_z_h));

            for j = 1:num_of_injections
                Bxh_noisy(:,:,j) = datah_b.Bx(:,j) + noise_time_series_x_h( ...
                    1:num_of_mdeit_channels_per_dim, :, j);

                Byh_noisy(:,:,j) = datah_b.By(:,j) + noise_time_series_y_h( ...
                    1:num_of_mdeit_channels_per_dim, :, j);

                Bzh_noisy(:,:,j) = datah_b.Bz(:,j) + noise_time_series_z_h( ...
                    1:num_of_mdeit_channels_per_dim, :, j);

                Bxi_noisy(:,:,j) = datai_b.Bx(:,j) + noise_time_series_x_i( ...
                    1:num_of_mdeit_channels_per_dim, :, j);

                Byi_noisy(:,:,j) = datai_b.By(:,j) + noise_time_series_y_i( ...
                    1:num_of_mdeit_channels_per_dim, :, j);

                Bzi_noisy(:,:,j) = datai_b.Bz(:,j) + noise_time_series_z_i( ...
                    1:num_of_mdeit_channels_per_dim, :, j);
            end

            % Now we make a decision on how to convert the time-series into
            % a value. For now just take the average

            Bxh_noisy = squeeze(mean(Bxh_noisy,2));
            Byh_noisy = squeeze(mean(Byh_noisy,2));
            Bzh_noisy = squeeze(mean(Bzh_noisy,2));

            Bxi_noisy = squeeze(mean(Bxi_noisy,2));
            Byi_noisy = squeeze(mean(Byi_noisy,2));
            Bzi_noisy = squeeze(mean(Bzi_noisy,2));
            


            % At this point, lets do the


            % DEBUG
            % Bxh_corrected = spherical_harmonic_based_noise_rejection( Bxh_noisy ,noise_structure);
            % Bxi_corrected = spherical_harmonic_based_noise_rejection( Bxi_noisy ,noise_structure);
            % 
            % 
            % Byh_corrected = spherical_harmonic_based_noise_rejection( Byh_noisy ,noise_structure);
            % Byi_corrected = spherical_harmonic_based_noise_rejection( Byi_noisy ,noise_structure);
            % 
            % Bzh_corrected = spherical_harmonic_based_noise_rejection( Bxh_noisy ,noise_structure);
            % Bzi_corrected = spherical_harmonic_based_noise_rejection( Bxi_noisy ,noise_structure);
            



            % % Add the noise at random time in time series
            % id = randi(size(noise_time_series_x,2));
            % Bxi_noisy = datai_b.Bx + squeeze(noise_time_series_x( ...
            %     1:num_of_mdeit_channels_per_dim, id, 1:num_of_injections));
            %
            % id = randi(size(noise_time_series_x,2));
            % Byi_noisy = datai_b.By + squeeze(noise_time_series_y( ...
            %     1:num_of_mdeit_channels_per_dim, id, 1:num_of_injections));
            %
            % id = randi(size(noise_time_series_x,2));
            % Bzi_noisy = datai_b.Bz + squeeze(noise_time_series_z( ...
            %     1:num_of_mdeit_channels_per_dim, id, 1:num_of_injections));
            

            dBx = datai_b.Bx(:)-datah_b.Bx(:);
            dBy = datai_b.By(:)-datah_b.By(:);
            dBz = datai_b.Bz(:)-datah_b.Bz(:);
                
            % FOR NOW, LETS ADD NOISE ONLY COMING FROM THE MEASUREMENT OF
            % Bzi, BUT THIS IS NOT VERY CORRECT!!!!!!!!!!!!!!!!!

            dBx_noisy = Bxi_noisy(:) - datah_b.Bx(:);
            dBy_noisy = Byi_noisy(:) - datah_b.By(:);
            dBz_noisy = Bzi_noisy(:) - datah_b.Bz(:);
            

            % Store the noise
            noise_b_x = dBx_noisy - dBx;
            noise_b_y = dBy_noisy - dBy;
            noise_b_z = dBz_noisy - dBz;


            % dBx_noisy = Bxi_noisy(:) - Bxh_noisy(:);
            % dBy_noisy = Byi_noisy(:) - Byh_noisy(:);
            % dBz_noisy = Bzi_noisy(:) - Bzh_noisy(:);

            dB = [dBx,dBy,dBz];
            dB_noisy = [dBx_noisy,dBy_noisy,dBz_noisy];

            snr = 0;
            for i = 1:3
                if not(all(dB(:,i)<1e-12))
                    snr = snr + 20*log(...
                        signal_amplitude(dB(:,i))/signal_amplitude(dB_noisy(:,i)-dB(:,i))...
                        );
                end
            end

            snr = snr/sum([any(dB(:,1)>1e-12),any(dB(:,2)>1e-12),any(dB(:,3)>1e-12)]);

            % Add white noise to EIT such that SNR is the same
            warning('Noise case: OPM. Adding white noise to EIT anyway');
            
            du = datai_u.meas-datah_u.meas;
            noise_structure.num_of_channels = num_of_eit_channels;
            noise_amplitude_eit = signal_amplitude(du) / 10^(snr/20);
            
            noise_u = randn(size(du,1),1);
            noise_u = noise_amplitude_eit*noise_u/signal_amplitude(noise_u);
            du_noisy = du + noise_u;
        end
        
        % Output
        out_data_b.dBx = dBx_noisy;
        out_data_b.dBy = dBy_noisy;
        out_data_b.dBz = dBz_noisy;

        out_data_u.meas_du = du_noisy;

        noise_b.noise_b_x = noise_b_x;
        noise_b.noise_b_y = noise_b_y;
        noise_b.noise_b_z = noise_b_z;
        

    case 'absolute'

        error('Not implemented yet')

end

end




%Article: Spherical harmonic based noise rejection and neuronal sampling with multi-axis OPMs
function B_corrected = spherical_harmonic_based_noise_rejection(B_field,noise_structure)
assert(isfield(noise_structure,'D'),'noise_structure must have a field "D"');

D = noise_structure.D;
% Spherical harmonic based noise rejection and neuronal sampling with multi-axis OPMs
S=[];
S.D = D;
S.scale = 1;
S.li = 11; % order of harmonics

S.reg= false;  % irregular spherical harmonics (neural space)
A = spm_opm_vslm(S);

S.reg = true; %regular spherical harmonics (interference space)
B = spm_opm_vslm(S);

% Get the channels which haven't been subjected to feedback (only 22!)
feedback_channel_ids = selectchannels(D, 'regexp_G2-[A|D].*');
feedback_channels = intersect(D.chanlabels(feedback_channel_ids),D.sensors('MEG').label);

meginds = intersect(indchannel(D, D.sensors('MEG').label), ...
    indchantype(D, 'MEGMAG', 'GOOD')); % all good meg sensors ids

feedbackinds = intersect(meginds,indchannel(D, feedback_channels));

meginds = setdiff(meginds, feedbackinds);

all_meg_ids = indchannel(D, D.sensors('MEG').label);
[~, ids] = ismember(meginds, all_meg_ids); %ids gives the index in all_meg_ids where each meginds is found.

A = A(ids,:); %select only non feedback channels

B = B(ids,:); %select only non feedback channels

% Number of lines of these matrices is numel(sensors(D,'MEG')). This is
% using the feedback channels as well!

%Concatenate columnwise
T = [A,B];

%THERE ARE SEVERAL SOLUTIONS!
c = pinv(T)*B_field; 

% c2 = T\B_field;
% c3 = pcg(T'*T+0*eye(size(T'*T)),T'*B_field(:,1));
% c4 = linsolve(T,B_field(:,1));

% Project data into subspace of A basis
c_a = c(1:size(c,1)/2,:);
c_b = c(size(c,1)/2+1:end,:);

B_corrected = A*c_a;

end


% Generate noise time-series to add to generate data
function [noise,t_id] = generate_noise(noise_structure)

if nargin<3
    N = [];
end

switch lower(noise_structure.type)
    case 'white'
        noise =  generate_white_noise(noise_structure);
        t_id = [];
    case 'pink'
        noise = generate_pink_noise(noise_structure);
        t_id = [];
    case 'brown'
        noise = generate_brown_noise(noise_structure);
        t_id = [];
    case 'opm'
        [noise,t_id] = generate_opm_noise(noise_structure,t_id);
    otherwise
        error('Unknown noise type: %s. Use "white", "pink", "brown" or "opm".', type);
end
end

%% generate_opm_noise
function [noise,t_id] = generate_opm_noise(noise_structure,t_id_0)

assert(isfield(noise_structure,'measurement_protocol'));

if not(strcmp(noise_structure.measurement_protocol,'default_rui'))
    error('Need to implement another measurement protocol!')
end

% This is data from a particular noise sampling from article: Data from Real-time, model-based magnetic field correction for moving, wearable MEG
D = noise_structure.D;

time = D.time;

assert(isfield(noise_structure,'epoch_time'),'noise_structure must have field "epoch_time"');
epoch_time = noise_structure.epoch_time;

assert(isfield(noise_structure,'num_of_measurements'),'noise_structure must have field "num_of_measurements"');
num_of_measurements = noise_structure.num_of_measurements;

assert(isfield(noise_structure,'B0'),'noise_structure must have field "H0"');
B0 = noise_structure.B0;


assert(isfield(noise_structure,'num_of_channels'),'noise_structure must have field "num_of_channels"');
num_of_channels = noise_structure.num_of_channels;

if noise_structure.num_of_channels~=22
    error('OPM noise generation only working for 22 channels for now! Need to research synthetic noise generation for OPMs if we allow more channels than this')
end

% Get the channels which haven't been subjected to feedback (only 22!)
feedback_channel_ids = selectchannels(D, 'regexp_G2-[A|D].*');
feedback_channels = intersect(D.chanlabels(feedback_channel_ids),D.sensors('MEG').label);

meginds = intersect(indchannel(D, D.sensors('MEG').label), ...
    indchantype(D, 'MEGMAG', 'GOOD')); % all good meg sensors ids

feedbackinds = intersect(meginds,indchannel(D, feedback_channels));

meginds = setdiff(meginds, feedbackinds);

% Map each channel in 1:num_of_channels to a random index in
% meginds

% Randomly select indices from meginds for each element of
% 1:num_of_channels

% Save current RNG state
state = rng;

seed = sum(double(meginds)) + num_of_channels * 100000;
rng(seed);
rand_ids = randi(length(meginds), length(1:num_of_channels), 1);
% Map each element of A to a value in B
mapped_ids = meginds(rand_ids);

% Restore original RNG state
rng(state);


% % Spherical harmonic based noise rejection and neuronal sampling with multi-axis OPMs
% S=[];
% S.D = D;
%
% S.scale = 1e-3;
% S.li = 11; % order of harmonics
%
% S.reg= false;  % irregular spherical harmonics (neural space)
% A = spm_opm_vslm(S);
%
% S.reg = true; %regular spherical harmonics (interference space)
% B = spm_opm_vslm(S);
%
% all_meg_ids = indchannel(D, D.sensors('MEG').label);
% [~, ids] = ismember(meginds, all_meg_ids); %ids gives the index in all_meg_ids where each meginds is found.
%
% A = A(ids,:); %select only non feedback channels
%
% B = B(ids,:); %select only non feedback channels
%
% % Number of lines of these matrices is numel(sensors(D,'MEG')). This is
% % using the feedback channels as well!
%
% %Concatenate columnwise
% T = [A,B];

% Pick a random epoch of time for each channel
dt = time(2) - time(1);          % sampling interval
n_samples_epoch = round(epoch_time / dt); % number of samples in the epoch

% Pick starting time at random, for each channel/measurement such that the epoch fits the
% time window
% Generate random integers from n_samples_epoch+1 to nsamples(D)-n_samples_epoch
% start_time_ids = ...
%     n_samples_epoch+randi(nsamples(D)-2*n_samples_epoch,num_of_channels,num_of_measurements);
% if any(start_time_ids(:)>2401480) || any(start_time_ids(:)<1)
%     error('This should not happen!');
% end

noise = zeros(num_of_channels,n_samples_epoch,num_of_measurements);

% Sample noise according to the following measurement protocol:
% For each injection pattern, measure the time-series of all the OPMS
% simultaneously on a window of time called the epoch. Then do the same for
% the next injection pattern for the next epoch time window.

%The start index t_id should be random!
t_id_start = t_id_0;

% figure
% hold on
% Dplot = D(mapped_ids,:,1)';
% plot(Dplot);
% plot(ones(100,1)*t_id_0,linspace(min(Dplot(:)),max(Dplot(:))),'r-')
% hold off

for j = 1:num_of_measurements
    t_id = t_id_start + (j-1)*n_samples_epoch+1;
    % Sample noise time series from data for each channel
    N = D(mapped_ids,t_id:t_id+n_samples_epoch-1,1);

    noise(:,:,j) = N;
end

t_id = t_id_start +  num_of_measurements*n_samples_epoch;



%     for i = 1:num_of_channels
% for j = 1:num_of_measurements
%         % Sample noise time series from data for each channel ( the time
%         % window also depends on the channel, so we need a smarter way if
%         % we want to vectorize)
%         N = D(mapped_ids(i),start_time_ids(i,j):start_time_ids(i,j)+n_samples_epoch-1,1);
%
%         noise(i,:,j) = N;
%     end
% end



% Indeed if we do SSS for the noise, N_a is bigger than N, but N_a + N_b =
% N up to some numerical error !!! So this is working properly, but must be
% applied for the generated magnetic fields probably, can't be applied to
% noise.
%
% for j = 1:num_of_measurements
%
%     N = D(meginds,start_time_ids(1,j):start_time_ids(1,j)+n_samples_epoch-1,1);
%
%     % Solve for coefficients in irreg/reg spherical harmonic basis
%     % (Moore-Penrose pseudo-inverse seems to be what works best.
%     % Regularized pcg solution is the same when lambda is 1e-7)
%     c = pinv(T)*N;
%
%     % Project data into subspace of A basis
%     c_a = c(1:size(c,1)/2,:);
%     % c_b = c(size(c,1)/2+1:end,:);
%
%
%     N_a = A*c_a;
%     % N_b = B*c_b;
% end



% Rescale noise to adimensional units ( result comes in fT)
noise = noise*1e-15/B0;

end

%% generate_pink_noise
function noise = generate_pink_noise(noise_structure)

assert(isfield(noise_structure,'epoch_time'),'noise_structure must have field "epoch_time"');
n_samples_epoch = noise_structure.n_samples_epoch;

assert(isfield(noise_structure,'num_of_channels'),'noise_structure must have field "num_of_channels"');
num_of_channels = noise_structure.num_of_channels;


%Power spectral density ~ 1/f

% Generate white Gaussian noise
w = randn(num_of_channels,n_samples_epoch);

% FFT
W =  fft(w,[],2);

% Frequency vector
freqs = (0:n_samples_epoch-1)'/n_samples_epoch;  % normalized frequency

% Avoid division by zero at DC
freqs(1) = freqs(2);

% 1/sqrt(f) filter magnitude (since PSD ∝ |H(f)|^2)
H = 1 ./ sqrt(freqs);

% Apply coloring filter
X = W .* repmat(H.', size(W,1), 1);

% Inverse FFT → pink noise
x = real(ifft(X,[],1));

% Normalize noise amplitude
noise = x ./ signal_amplitude(x').';
end

%% generate_brown_noise
function noise = generate_brown_noise(noise_structure)

assert(isfield(noise_structure,'epoch_time'),'noise_structure must have field "epoch_time"');
n_samples_epoch = noise_structure.n_samples_epoch;

assert(isfield(noise_structure,'num_of_channels'),'noise_structure must have field "num_of_channels"');
num_of_channels = noise_structure.num_of_channels;

% Brown noise = cumulative sum of white noise
w = randn(num_of_channels,n_samples_epoch);
noise = cumsum(w,2);


% Normalize noise  
noise = noise ./ signal_amplitude(noise').';

end




%% generate_white_noise
function noise = generate_white_noise(noise_structure)

assert(isfield(noise_structure,'epoch_time'),'noise_structure must have field "epoch_time"');
n_samples_epoch = noise_structure.n_samples_epoch;

assert(isfield(noise_structure,'num_of_channels'),'noise_structure must have field "num_of_channels"');
num_of_channels = noise_structure.num_of_channels;

n = randn(num_of_channels,n_samples_epoch);

noise = n./signal_amplitude(n').';
end
%% validate_options
function [noise_structure,D] = validate_options(options)

assert(isfield(options,'noise_structure'),'options must have a "noise_structure" field');
noise_structure = options.noise_structure;

valid_types = {'white','pink','brown','opm'};
assert(any(strcmp(noise_structure.type,valid_types)),'not a valid type');

assert(isfield(noise_structure,'type'),'noise_structure must have a "type" field');
assert(isfield(noise_structure,'snr'),'noise_structure must have a "snr" field');


% Load necessary data if noise type is opm
if strcmp(noise_structure.type,'opm')
    assert(isfield(options,'data_path'),'options must have a "data_path" field');

    %Check if file exists
    if exist(options.data_path,'file')
        D = spm_eeg_load(options.data_path);
        noise_structure.D = D;
    else
        error('%s not found',options.data_path)
    end
end

end

%% This function is used to switch between definitios of signal amplitude
function amplitudes = signal_amplitude(s)

if size(s,1) == 1 && size(s,2) > 1 %if line vector, convert to column vector
    s = s(:);
end

% % % Define the amplitude as the std for each column
% amplitudes = std(s,1);

% Define the amplitude as largest deviation to mean
amplitudes = zeros(1,size(s,2));
means = mean(s,1);

for j = 1:size(s,2)
    amplitudes(j) = norm(s(:,j)-means(j),'inf');
end

end


function snr = compute_snr(signal,noise)
    snr = 20.*log10(signal_amplitude(signal)./signal_amplitude(noise));
end