%% For 3-axis MDEIT
function dB = computeSensitivities3axis(img,lambdatimesdAdp,lambdatimesdbdp,A,delta,eidorsFolder)

% fprintf('Computing Jacobian:\n')

numStim = size(img.fwd_model.stimulation,2);
numElements = size(img.fwd_model.elems,1);

sigma = img.elem_data;

%Get EIT forward solution
uinh = fwd_solve(img);
u = uinh.volt;

[GammaX,GammaY,GammaZ] = computeGammaMatrices(img,sigma);

numOfSensors = size(GammaX,1);

dGammaXcell = img.dGammaXcell;
dGammaYcell = img.dGammaYcell;
dGammaZcell = img.dGammaZcell;

dftempx= zeros(numStim,numOfSensors,numElements);
dftempy= zeros(numStim,numOfSensors,numElements);
dftempz= zeros(numStim,numOfSensors,numElements);

%USE PARFOR HERE For parallelization
parfor m = 1:numOfSensors

    % fprintf('Magnetometer %i\n',m);

    dfdx = @(x,p) GammaX(m,:)';
    dfdy = @(x,p) GammaY(m,:)';
    dfdz = @(x,p) GammaZ(m,:)';

    % dLmdu = @(x,p)

    for n = 1:numStim
        un = u(:,n);

        dBxdsigma = un'*dGammaXcell{m};
        dBydsigma = un'*dGammaYcell{m};
        dBzdsigma = un'*dGammaZcell{m};

        dfxdp = @(u,p) dBxdsigma;
        dfydp = @(u,p) dBydsigma;
        dfzdp = @(u,p) dBzdsigma;

        M = @(sigma) A; %THIS IS NOT COMPUTING THE MATRIX! IT COMES PRE-COMPUTED, BECAUSE OF A PROBLEM WITH PARALLELIZATION!!!!

        dftempx(n,m,:) =  adjoint(un,sigma,M,...
            lambdatimesdAdp,lambdatimesdbdp,...
            dfdx,dfxdp,delta);

        dftempy(n,m,:) = adjoint(un,sigma,M,...
            lambdatimesdAdp,lambdatimesdbdp,...
            dfdy,dfydp,delta);

        dftempz(n,m,:) =  adjoint(un,sigma,A,...
            lambdatimesdAdp,lambdatimesdbdp,...
            dfdz,dfzdp,delta);
    end
end

%Convert df to be a matrix
dB = zeros(3*numStim*numOfSensors,numElements);
for n = 1:numStim
    for m = 1:numOfSensors

        idx = (1-1)*numStim*numOfSensors + (n-1)*numOfSensors+m;
        idy = (2-1)*numStim*numOfSensors + (n-1)*numOfSensors+m;
        idz = (3-1)*numStim*numOfSensors + (n-1)*numOfSensors+m;

        dB(idx,:) = dftempx(n,m,:);
        dB(idy,:) = dftempy(n,m,:);
        dB(idz,:) = dftempz(n,m,:);
    end
end

end
