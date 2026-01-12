function lambda = LCurveMethod(x0,res,jac)

lambdaVec = logspace(-2,1,5);

residualNorms = zeros(1,length(lambdaVec));
xNorms = zeros(1,length(lambdaVec));

maxNumOfIterations = 1;
tol = 1e-5;

for i = 1:length(lambdaVec)

    x = x0;

    % 0th Tykhonov prior
    RtRprior = @(Jk) speye(length(x),length(x));

    Wsqrt = [];

    [x,~,~] = ...
        gaussNewtonNoisy(res,jac,x,tol,maxNumOfIterations,Wsqrt,RtRprior,lambdaVec(i));
    
    % [x,~,~] = ...
    %     gaussNewtonV3(res,jac,x,tol,maxNumOfIterations,Wsqrt,RtRprior,lambdaVec(i));
    residualNorms(i) = norm(res(x),2);
    xNorms(i) = norm(x,2);
end

%% Compute the curvature according to L-curve article



%%
figure
hold on
plot(residualNorms,xNorms,'bo-');

xlabel('$|| r(\sigma_{\lambda}) ||_2$','Interpreter','latex');
ylabel('$|| \sigma_{\lambda} ||_2$','Interpreter','latex');

step = max(round(length(lambdaVec)/3),1);
for i = step:step:length(lambdaVec)
    str = sprintf('l = %.2d ',lambdaVec(i));
    text(1.01*residualNorms(i),xNorms(i),str,'Interpreter','latex')
end

title('L-curve','Interpreter','latex')
box on
grid on;grid minor
set(gca,'YScale','log')
set(gca,'XScale','log')

end