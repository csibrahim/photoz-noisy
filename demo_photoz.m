
rng(1); % fix random seed

addpath GPz/ % path to GPz
addpath(genpath('minFunc_2012/'))       % path to minfunc

if(false) % set to false to set the options externally 
    
    rng(1); % fix random seed

    m = 10;

    method = 'VD';

    heteroscedastic = true;

    normalize = true;

    maxIter = 500;
    maxAttempts = 100;


    trainSplit = 0.2;
    validSplit = 0.2;
    testSplit  = 0.6;

    csl_method = 'normal';

    trainOption = 3;    % 1=Features, 2=Sampling, 3=Parameters and 4=Weights
    samplingSize = 30;  % if Sampling is used

    transform = 3;   % Assume 1 = log-normal, 2=gamma, 3=luptitudes and otherwise linear
     
    trainPath = 'data/XMM_data_PH_27_11_18.csv';
    testPath = 'data/COSMOS_data_PH_27_11_18.csv';
    outPath = [];
end

%%%%%%%%%%%%%% Prepare data %%%%%%%%%%%%%% 
                                        

X = csvread(trainPath);
Y = X(:,end);
ra = X(:,1);
dec = X(:,2);

X = X(:,3:26);
 
remove = sum(X==0,2)>0;

X(remove,:)=[];
Y(remove)=[];

[n,d] = size(X);
filters = d/2;
  
% get the weights for cost-sensitive learning
omega = getOmega(Y,csl_method,[]);

% select training, validation and testing sets from the data
[training,validation,testing] = sample(n,trainSplit,validSplit,testSplit); 

Psi = X(:,filters+1:end).^2;
X(:,filters+1:end) = [];

switch(transform)
    case 1 % log-normal
        Psi = log(X.^2+Psi)-2*log(X);
        X = log(X)-0.5*Psi;

    case 2 % gamma

        theta = Psi./X;
        k = X./theta;
        X = psi(k)+log(theta);
        Psi = psi(1,k);

    case 3 % luptitudes

        a = 2.5/log(10);
        b = mode(sqrt(Psi));% could also be mean or median
        Xn = bsxfun(@rdivide,X,2*b);

        Psi = (a)^2*bsxfun(@rdivide,Psi./(1+Xn.^2),4*b.^2);
        X = -a*bsxfun(@plus,asinh(Xn),log(b));
end

if(trainOption == 1) % use errors as features
    
    X = [X Psi];
    
    Psi = [];
    
elseif(trainOption==2)% use erorrs to sample
    
    X = reshape(repmat(X',samplingSize,1),filters,samplingSize*n)';
    Psi = reshape(repmat(Psi',samplingSize,1),filters,samplingSize*n)';
    
    X = randn(samplingSize*n,filters).*sqrt(Psi)+X;
    Y = reshape(repmat(Y',samplingSize,1),samplingSize*n,1);
    omega = reshape(repmat(omega',samplingSize,1),samplingSize*n,1);
    
    training = reshape(repmat(training',samplingSize,1),samplingSize*n,1);
    validation = reshape(repmat(validation',samplingSize,1),samplingSize*n,1);
    testing = reshape(repmat(testing',samplingSize,1),samplingSize*n,1);
    
    Psi = [];
elseif(trainOption==3) % treat the mag-errors as input noise variance
    % if parameters do nothing
elseif(trainOption == 4) % use errors as weights in CSL
    
    sumLogPsi = -0.5*sum(log(Psi),2);
    
    omega = omega.*exp(sumLogPsi-max(sumLogPsi));
    
    Psi = [];
    
end



%%%%%%%%%%%%%% Fit the model %%%%%%%%%%%%%%

% initialize the model
model = init(X,Y,method,m,'omega',omega,'training',training,'heteroscedastic',heteroscedastic,'normalize',normalize,'Psi',Psi);
% train the model
model = train(model,X,Y,'omega',omega,'training',training,'validation',validation,'maxIter',maxIter,'maxAttempts',maxAttempts,'Psi',Psi); 


%%%%%%%%%%%%%% Compute Metrics %%%%%%%%%%%%%%

X = csvread(testPath);
Y = X(:,end);
ra = X(:,1);
dec = X(:,2);

X = X(:,3:26);
 
remove = sum(X==0,2)>0;

X(remove,:)=[];
Y(remove)=[];

Psi = X(:,filters+1:end).^2;
X(:,filters+1:end) = [];

switch(transform)
    case 1 % log-normal
        Psi = log(X.^2+Psi)-2*log(X);
        X = log(X)-0.5*Psi;

    case 2 % gamma

        theta = Psi./X;
        k = X./theta;
        X = psi(k)+log(theta);
        Psi = psi(1,k);

    case 3 % luptitudes

        a = 2.5/log(10);
        b = mode(sqrt(Psi));% could also be mean or median
        Xn = bsxfun(@rdivide,X,2*b);

        Psi = (a)^2*bsxfun(@rdivide,Psi./(1+Xn.^2),4*b.^2);
        X = -a*bsxfun(@plus,asinh(Xn),log(b));
end

sumLogPsi = -0.5*sum(log(Psi),2);
detPsi = exp(sumLogPsi);

if(trainOption == 1) % use errors as features
    
    X = [X Psi];
    
    Psi = [];
    
elseif(trainOption==2)% use erorrs to sample
    
    X = reshape(repmat(X',samplingSize,1),filters,samplingSize*n)';
    Psi = reshape(repmat(Psi',samplingSize,1),filters,samplingSize*n)';
    
    X = randn(samplingSize*n,filters).*sqrt(Psi)+X;
    Y = reshape(repmat(Y',samplingSize,1),samplingSize*n,1);

    Psi = [];
elseif(trainOption==3) % treat the mag-errors as input noise variance
    % if parameters do nothing
elseif(trainOption == 4) % use errors as weights in CSL
    
    Psi = [];
    
end

% use the model to generate predictions for the test set
[mu,sigma,nu,beta_i,gamma] = predict(X,model,'Psi',Psi);

n = length(mu);

if(trainOption==2)
    mu = reshape(mu,samplingSize,n/samplingSize);
    
    nu = reshape(nu,samplingSize,n/samplingSize);
    beta_i = reshape(beta_i,samplingSize,n/samplingSize);
    gamma = reshape(gamma,samplingSize,n/samplingSize);
    
    Y = reshape(Y,samplingSize,length(Y)/samplingSize)';
    Y = Y(:,1);
    
    mu2 = mean(mu.^2)';
    mu = mean(mu)';
    varMu = mu2-mu.^2;
    
    nu = mean(nu)';
    beta_i = mean(beta_i)';
    gamma = mean(gamma)'+varMu;
    
    sigma = nu+beta_i+gamma;
end

% mu     = the best point estimate
% nu     = variance due to data density
% beta_i = variance due to output noise
% gamma  = variance due to input noise
% sigma  = nu+beta_i+gamma

 
% save output as comma separated values
if(~isempty(outPath))
    csvwrite(outPath,[Y mu sigma nu beta_i gamma]);
end
