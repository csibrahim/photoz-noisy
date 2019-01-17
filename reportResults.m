
if(true) % set to false to set it externally
    folder = 'results/';
end

%%%%% preparing the variables %%%%%
files = dir([folder,'*.csv']);
names = {files.name};

k = numel(names);

data = csvread([folder,names{1}]);

n = size(data,1);

mus = zeros(n,k);
sigmas  = zeros(n,k);
nus = zeros(n,k);
beta_is = zeros(n,k);
gammas = zeros(n,k);

rmse = zeros(n,k);
mll  = zeros(n,k);
fr15 = zeros(n,k);
fr05 = zeros(n,k);
bias = zeros(n,k);

%%%%% compute the metrics and plot the scatter plot for each file %%%%%%
for i=1:k
    
    data = csvread([folder,names{i}]);
    Y = data(:,1); mus(:,i) = data(:,2); sigmas(:,i) = data(:,3); nus(:,i) = data(:,4); beta_is(:,i) = data(:,5); gammas(:,i) = data(:,6);
    
    %root mean squared error, i.e. sqrt(mean(errors^2))
    rmse(:,i) = sqrt(metrics(Y,mus(:,i),sigmas(:,i),@(y,mu,sigma) (y-mu).^2));

    % mean log likelihood mean(-0.5*errors^2/sigma -0.5*log(sigma)-0.5*log(2*pi))
    mll(:,i) = metrics(Y,mus(:,i),sigmas(:,i),@(y,mu,sigma) -0.5*(y-mu).^2./sigma - 0.5*log(sigma)-0.5*log(2*pi));

    % fraction of data where |z_spec-z_phot|/(1+z_spec)<0.15
    fr15(:,i) = metrics(Y,mus(:,i),sigmas(:,i),@(y,mu,sigma) 100*(abs(y-mu)./(y+1)<0.15));

    % fraction of data where |z_spec-z_phot|/(1+z_spec)<0.05
    fr05(:,i) = metrics(Y,mus(:,i),sigmas(:,i),@(y,mu,sigma) 100*(abs(y-mu)./(y+1)<0.05));

    % bias, i.e. mean(errors)
    bias(:,i) = metrics(Y,mus(:,i),sigmas(:,i),@(y,mu,sigma) y-mu);

    % now plot the density and the uncertainty plots side by side
%     figure;
    
%     set(gcf,'NumberTitle','off')
%     set(gcf,'Name',names{i})

    subplot(k,2,2*(i-1)+1);
    heat(Y,mus(:,i),log(sigmas(:,i)),0.01,0.1);title('Uncertainty');
    xlabel('Spectroscopic Redshift');
    ylabel('Photometric Redshift');
    colormap jet;
    
    subplot(k,2,2*(i-1)+2);
    heat(Y,mus(:,i),[],0.01,0.1);
    title('Density');
    xlabel('Spectroscopic Redshift');
    ylabel('Photometric Redshift');
    colormap jet;

end

pos = get(gcf,'Position');
set(gcf,'Position',[pos(1)/4  pos(2)   2*pos(3)   2*pos(4)]);
drawnow;


% print metrics for the entire data for each file
fprintf('RMSE\t\tMLL\t\tFR15\t\tFR05\t\tBIAS\n')
fprintf('%f\t%f\t%f\t%f\t%f\n',[rmse(end,:);mll(end,:);fr15(end,:);fr05(end,:);bias(end,:)])

% plot the change in metrics as functions of data percentage for each file
x = [1 5:5:100];
ind = round(x*n/100);

figure;plot(x,rmse(ind,:),'o-');xlabel('Percentage of Data');ylabel('RMSE');legend(names);drawnow
figure;plot(x,mll(ind,:),'o-');xlabel('Percentage of Data');ylabel('MLL');legend(names);drawnow
figure;plot(x,fr05(ind,:),'o-');xlabel('Percentage of Data');ylabel('FR05');legend(names);drawnow
figure;plot(x,fr15(ind,:),'o-');xlabel('Percentage of Data');ylabel('FR15');legend(names);drawnow
figure;plot(x,bias(ind,:),'o-');xlabel('Percentage of Data');ylabel('BIAS');legend(names);drawnow

% plot the mean and the standard deviation of different scores as functions of spectroscopic redshift using 20 bins
[centers,means,stds] = bin(Y,-bsxfun(@minus,mus,Y),20);
figure;errorbar(repmat(centers,1,k),means,stds,':','LineWidth',2);xlabel('Spectroscopic Redshift');ylabel('Bias');legend(names);drawnow

[centers,means,stds] = bin(Y,sqrt(nus),20);
figure;errorbar(repmat(centers,1,k),means,stds,':','LineWidth',2);xlabel('Spectroscopic Redshift');ylabel('Model Uncertainty');legend(names);drawnow

[centers,means,stds] = bin(Y,sqrt(beta_is),20);
figure;errorbar(repmat(centers,1,k),means,stds,':','LineWidth',2);xlabel('Spectroscopic Redshift');ylabel('Noise Uncertainty');legend(names);drawnow

[centers,means,stds] = bin(Y,sqrt(gammas),20);
figure;errorbar(repmat(centers,1,k),means,stds,':','LineWidth',2);xlabel('Spectroscopic Redshift');ylabel('Input Noise Uncertainty');legend(names);drawnow
