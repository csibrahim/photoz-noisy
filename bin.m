function [centers,means,stds] = bin(x,y,bins,xmax)

    if(nargin==2)
        bins = 100;
    end
    if(nargin==3)
        xmax = max(x)
    end

    centers = linspace(0,xmax,bins)';
    n = length(x);
    m = length(centers);
    [~,id] = min(Dxy(x,centers),[],2);
    counts = full(sum(sparse(1:length(x),id,1,n,m)))';

    remove = counts==0;
    counts(remove) = 1;

    k = size(y,2);
    means = zeros(length(centers),k);
    stds = zeros(length(centers),k);

    for i=1:k
        sums = full(sum(sparse(1:length(x),id,y(:,i),n,m)))';

        means(:,i) = sums./counts;

        sumssqrs = full(sum(sparse(1:length(x),id,y(:,i).^2,n,m)))';

        stds(:,i) = sumssqrs./counts;
        stds(:,i) = sqrt(stds(:,i)-means(:,i).^2);

    end
    
    means(remove,:) = [];
    stds(remove,:) = [];
    centers(remove,:) = [];
    
end

