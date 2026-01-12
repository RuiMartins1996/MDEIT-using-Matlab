function dB = ...
    computeSensitivitiesNoParallel(img,lambdatimesdAdp,lambdatimesdbdp,A,delta,eidorsFolder)

fprintf('Computing Jacobian (no parallelization):\n')

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

for m = 1:numOfSensors

    fprintf('Magnetometer %i\n',m);

    dfdx = @(x,p) GammaX(m,:)';
    dfdy = @(x,p) GammaY(m,:)';
    dfdz = @(x,p) GammaZ(m,:)';

    for n = 1:numStim
        un = u(:,n);

        dBxdsigma = un'*dGammaXcell{m};
        % dBydsigma = un'*dGammaYcell{m};
        % dBzdsigma = un'*dGammaZcell{m};

        dfdp = @(u,p) dBxdsigma;
        % dfdpy = @(u,p) dBydsigma;
        % dfdpz = @(u,p) dBzdsigma;
        
        M = @(sigma) A; %THIS IS NOT COMPUTING THE MATRIX! IT COMES PRE-COMPUTED, BECAUSE OF A PROBLEM WITH PARALLELIZATION!!!!

        dftempx(n,m,:) =  adjoint(un,sigma,M,...
            lambdatimesdAdp,lambdatimesdbdp,...
            dfdx,dfdp,delta);

        % dftempy(n,m,:) =  adjoint(un,sigma,A,...
        %     lambdatimesdAdp,lambdatimesdbdp,...
        %     dfdy,dfdpy,delta);
        % 
        % dftempz(n,m,:) =  adjoint(un,sigma,A,...
        %     lambdatimesdAdp,lambdatimesdbdp,...
        %     dfdz,dfdpz,delta);
    end
end

%Convert df to be a matrix
dB = zeros(numStim*numOfSensors,numElements);
for n = 1:numStim
    for m = 1:numOfSensors

            id = numOfSensors*(n-1)+(m-1)+1;

            dB(id,:) = dftempx(n,m,:);
    end
end
end

