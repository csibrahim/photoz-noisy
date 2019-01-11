function [centers,means,stds] = bin_single(x,y,bins,xmax)

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

    sums = full(sum(sparse(1:length(x),id,y,n,m)))';
    means = sums./counts;

    sumssqrs = full(sum(sparse(1:length(x),id,y.^2,n,m)))';
    stds = sumssqrs./counts;
    stds = sqrt(stds-means.^2);

    means(remove) = [];
    stds(remove) = [];
    centers(remove) = [];
end

