rng(1); % fix random seed
addpath GPz/ % path to GPz
addpath(genpath('minFunc_2012/'))       % path to minfunc
 
m = 100;
 
method = 'VD';
 
heteroscedastic = true;
 
normalize = true;

maxIter = 500;
maxAttempts = 100;
 
 
trainSplit = 0.2;
validSplit = 0.2;
testSplit  = 0.6;

csl_method = 'normal';

trainOption = 3;    % 1=Features, 2=Sampling and 3=Parameters
samplingSize = 30;  % if Sampling is used

transform = 0;   % assum a log-normal distribution and transform the mean and variance accordingly

binWidth = 0.1;

%%%%%%%%%%%%%% Prepare data %%%%%%%%%%%%%% 
 
% dataPath = 'data/COSMOS_data_PH_27_11_18.csv'; 
dataPath = 'data/XMM_data_PH_27_11_18.csv';
outPath = [];
                                        

X = csvread(dataPath);
Y = X(:,end);
ra = X(:,1);
dec = X(:,2);

X = X(:,3:26);
 
remove = sum(X==0,2)>0;

X(remove,:)=[];
Y(remove)=[];

[n,d] = size(X);
filters = d/2;
 
% you can also select the size of each sample
% [training,validation,testing] = sample(n,10000,10000,10000);
  
% get the weights for cost-sensitive learning
omega = getOmega(Y,csl_method,binWidth);

% select training, validation and testing sets from the data
[training,validation,testing] = sample(n,trainSplit,validSplit,testSplit); 

if(trainOption == 1) % use errors as features
    
    if(transform)
        X = log(X);
    end
    
    Psi = [];
    
elseif(trainOption==2)% use erorrs to sample
    
    Psi = X(:,filters+1:end).^2;
    X(:,filters+1:end) = [];
    
    if(transform)
        Psi = log(X.^2+Psi)-2*log(X);
        X = log(X)-0.5*Psi;
    end
    
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
            b = sqrt(Psi);
            Psi = (((2.5/log(10))./(2*b))).*(b./(1 + (X./(2*b))));
            X = -2.5/log(10) * (asinh(X./(2*b)) + log(b));
    end
else % use errors as weights in CSL
    
    Psi = X(:,filters+1:end).^2;
    X(:,filters+1:end) = [];
    
    if(transform)
        Psi = log(X.^2+Psi)-2*log(X);
        X = log(X)-0.5*Psi;
    end
    
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

% use the model to generate predictions for the test set
[mu,sigma,nu,beta_i,gamma] = predict(X,model,'Psi',Psi,'selection',testing);

if(trainOption==2)
    mu = reshape(mu,samplingSize,sum(testing)/samplingSize);
    sigma = reshape(sigma,samplingSize,sum(testing)/samplingSize);
    nu = reshape(nu,samplingSize,sum(testing)/samplingSize);
    beta_i = reshape(beta_i,samplingSize,sum(testing)/samplingSize);
    gamma = reshape(gamma,samplingSize,sum(testing)/samplingSize);
    
    Y = reshape(Y,samplingSize,length(Y)/samplingSize)';
    Y = Y(:,1);
    
    testing = reshape(testing,samplingSize,length(testing)/samplingSize)';
    testing = testing(:,1);
    
    mu2 = mean(mu.^2)';
    mu = mean(mu)';
    varMu = mu2-mu.^2;
    
    sigma = mean(sigma)'+varMu;
    nu = mean(nu)'+varMu;
    beta_i = mean(beta_i)'+varMu;
    gamma = mean(gamma)'+varMu;
end


% mu     = the best point estimate
% nu     = variance due to data density
% beta_i = variance due to output noise
% gamma  = variance due to input noise
% sigma  = nu+beta_i+gamma

% compute metrics 
 
%root mean squared error, i.e. sqrt(mean(errors^2))
rmse = sqrt(metrics(Y(testing),mu,sigma,@(y,mu,sigma) (y-mu).^2)); 
 
% mean log likelihood mean(-0.5*errors^2/sigma -0.5*log(sigma)-0.5*log(2*pi))
mll = metrics(Y(testing),mu,sigma,@(y,mu,sigma) -0.5*(y-mu).^2./sigma - 0.5*log(sigma)-0.5*log(2*pi));
 
% fraction of data where |z_spec-z_phot|/(1+z_spec)<0.15
fr15 = metrics(Y(testing),mu,sigma,@(y,mu,sigma) 100*(abs(y-mu)./(y+1)<0.15));
 
% fraction of data where |z_spec-z_phot|/(1+z_spec)<0.05
fr05 = metrics(Y(testing),mu,sigma,@(y,mu,sigma) 100*(abs(y-mu)./(y+1)<0.05));
 
% bias, i.e. mean(errors)
bias = metrics(Y(testing),mu,sigma,@(y,mu,sigma) y-mu);
 
% print metrics for the entire data
fprintf('RMSE\t\tMLL\t\tFR15\t\tFR05\t\tBIAS\n')
fprintf('%f\t%f\t%f\t%f\t%f\n',rmse(end),mll(end),fr15(end),fr05(end),bias(end))
 
%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%% 
 
% reduce the sample for efficient plotting
[x,y,color,counts]=reduce(Y(testing),mu,sigma,200);
 
figure;heat(Y(testing),mu,log(sigma),0.01,0.1);title('Uncertainty');xlabel('Spectroscopic Redshift');ylabel('Photometric Redshift');colormap jet;
figure;heat(Y(testing),mu,[],0.01,0.1);title('Density');xlabel('Spectroscopic Redshift');ylabel('Photometric Redshift');colormap jet;
 
% plot the change in metrics as functions of data percentage
x = [1 5:5:100];
ind = round(x*length(rmse)/100);
 
figure;plot(x,rmse(ind),'o-');xlabel('Percentage of Data');ylabel('RMSE');
figure;plot(x,mll(ind),'o-');xlabel('Percentage of Data');ylabel('MLL');
figure;plot(x,fr05(ind),'o-');xlabel('Percentage of Data');ylabel('FR05');
figure;plot(x,fr15(ind),'o-');xlabel('Percentage of Data');ylabel('FR15');
figure;plot(x,bias(ind),'o-');xlabel('Percentage of Data');ylabel('BIAS');
 
% plot mean and standard deviation of different scores as functions of spectroscopic redshift using 20 bins
[centers,means,stds] = bin(Y(testing),Y(testing)-mu,20);
figure;errorbar(centers,means,stds,'s');xlabel('Spectroscopic Redshift');ylabel('Bias');
 
[centers,means,stds] = bin(Y(testing),sqrt(nu),20);
figure;errorbar(centers,means,stds,'s');xlabel('Spectroscopic Redshift');ylabel('Model Uncertainty');
 
[centers,means,stds] = bin(Y(testing),sqrt(beta_i),20);
figure;errorbar(centers,means,stds,'s');xlabel('Spectroscopic Redshift');ylabel('Noise Uncertainty');

[centers,means,stds] = bin(Y(testing),sqrt(gamma),20);
figure;errorbar(centers,means,stds,'s');xlabel('Spectroscopic Redshift');ylabel('Input Noise Uncertainty');
 
% save output as comma separated values
if(~isempty(outPath))
    csvwrite([method,'_',num2str(m),'_',csl_method,'.csv'],[Y(testing) mu sigma nu beta_i gamma]);
end
