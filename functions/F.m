function [Bx,By,Bz] = F(img,sensorLocations)
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

%CAREFUL HERE: img.fwd_solve.get_all_meas = 1 must be set for the .volt
%field to exist

Bx = GammaX*u.volt;
By = GammaY*u.volt;
Bz = GammaZ*u.volt;

end

