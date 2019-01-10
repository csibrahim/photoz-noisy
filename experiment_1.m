
%%%%% set shared parameters %%%%%%

m = 100;
 
method = 'VD';
 
heteroscedastic = true;
 
normalize = true;

maxIter = 500;
maxAttempts = 50;
 
 
trainSplit = 0.8;
validSplit = 0.2;
testSplit  = 0;

csl_method = 'normal';

trainOption = 1;    % 1=Features, 2=Sampling, 3=Parameters and 4=Weights
samplingSize = 10;  % if sampling (trainOption=2) is used

transform = 3;   % Assume 1 = log-normal, 2=gamma, 3=luptitudes and otherwise linear

trainPath = 'data/XMM_data_PH_27_11_18.csv'; % data to build the model
testPath = 'data/COSMOS_data_PH_27_11_18.csv'; % data to test the model

folder = 'results/'; % folder to store the results


%%%%%% Change the desired variables and run the experiments %%%%%

trainOption = 1; % use noise as features
outPath = [folder,'trainOption_',num2str(trainOption),'.csv']; % where to store the results
demo_photoz; % run the experiemnt

trainOption = 2;% now use the sampling method
outPath = [folder,'trainOption_',num2str(trainOption),'.csv'];% where to store the results
demo_photoz;% run the experiemnt again with the changed parameter

trainOption = 3;% now use the parameter method
outPath = [folder,'trainOption_',num2str(trainOption),'.csv'];% where to store the results
demo_photoz;% run the experiemnt again with the changed parameter

trainOption = 4;% now use the noise as weights
outPath = [folder,'trainOption_',num2str(trainOption),'.csv'];% where to store the results
demo_photoz;% run the experiemnt again with the changed parameter

%%%%%%%% now plot the results %%%%%%%

reportResults

