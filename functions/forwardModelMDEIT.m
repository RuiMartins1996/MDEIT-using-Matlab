function [Bx,By,Bz] = forwardModelMDEIT(img,sensorLocations)

if nargin<2
    if ~isfield(img,'sensorLocations')
        error('Field sensorLocations must be present in img');
    end
    sensorLocations = img.sensorLocations;
end

numSensors = size(sensorLocations,1);
numStim = length(img.fwd_model.stimulation);

Bx = zeros(numSensors*numStim,1);
By = zeros(numSensors*numStim,1);
Bz = zeros(numSensors*numStim,1);

%% Forward solve the EIT model
u = fwd_solve( img );

%% FUNCTION: Compute the gamma matrices

    function [GammaX,GammaY,GammaZ] = computeGammaMatrices(img,sigma)

        mu0 = img.mu0;
        Omega = sparse(1:length(img.elem_volume),1:length(img.elem_volume),img.elem_volume);
        Sigma = sparse(1:length(sigma),1:length(sigma),sigma);
     
        RxTimesOmega = img.R.Rx*Omega;
        RyTimesOmega = img.R.Ry*Omega;
        RzTimesOmega = img.R.Rz*Omega;

        GammaX = -mu0/(4*pi)*RzTimesOmega*Sigma*img.G.Gy+...
            mu0/(4*pi)*RyTimesOmega*Sigma*img.G.Gz;

        GammaY = -mu0/(4*pi)*RxTimesOmega*Sigma*img.G.Gz+...
            mu0/(4*pi)*RzTimesOmega*Sigma*img.G.Gx;

        GammaZ = -mu0/(4*pi)*RyTimesOmega*Sigma*img.G.Gx+...
            mu0/(4*pi)*RxTimesOmega*Sigma*img.G.Gy;

    end

%% Output

[GammaX,GammaY,GammaZ] = computeGammaMatrices(img,img.elem_data);

for i = 1:numStim

    ui = u.volt(:,i);

    Bxi = GammaX*ui;
    Byi = GammaY*ui;
    Bzi = GammaZ*ui;

    idstart = (i-1)*numSensors+1;
    idend = (i-1)*numSensors+numSensors;

    Bx(idstart:idend) = Bxi;
    By(idstart:idend) = Byi;
    Bz(idstart:idend) = Bzi;
end
end

