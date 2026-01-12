function df = gradientAdjoint(sigma,ctx,options,Mfunc,lambdatimesdAdp,lambdatimesdbdp)

img = ctx.img;
Bxreal = ctx.Bxreal;

tol = options.tol;

dGammaXcell = img.dGammaXcell;
dGammaYcell = img.dGammaYcell;
dGammaZcell = img.dGammaZcell;

%Problem parameters
numElements = size(img.fwd_model.elems,1);
numStim = size(img.fwd_model.stimulation,2);
numSensors = length(dGammaXcell);

% Solve forward EIT problem
u = fwd_solve(img);

% Solve forward MDEIT problem
[GammaX,GammaY,GammaZ] = computeGammaMatrices(img,sigma);

Bx = GammaX*u.volt;
By = GammaY*u.volt;
Bz = GammaZ*u.volt;

% Get cost function gradient with adjoint method
df = zeros(1,numElements);

for i = 1:numStim
    ui = u.volt(:,i);

    ri = Bx(:,i)-Bxreal(:,i);
    
    dfdx = @(x,p) 2*ri'*GammaX;

    dGammaXTotal = zeros(numSensors,numElements);

    for m = 1:numSensors
        dGammaXTotal(m,:) =  ui'*dGammaXcell{m};
    end

    dfdsigma = 2*ri'*dGammaXTotal; 

    dfdp = @(x,p) dfdsigma;

    df = df + adjoint(ui,sigma,Mfunc,lambdatimesdAdp,lambdatimesdbdp,dfdx,dfdp,tol);
end

%Return a column vector
if size(df,2)>size(df,1)
    df = df';
    return
end


end

