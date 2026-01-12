%FUNCTION: Compute the gamma matrices. THIS SHOULD MAYBE BE STORED IN
%MEMORY AFTER THE FORWARD SOLVE, TO NOT REPEAT THE CALCULATION
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
