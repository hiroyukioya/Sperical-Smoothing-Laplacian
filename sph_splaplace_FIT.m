function [S,L, lambda] = sph_splaplace_FIT(Av, K,LapK, T, Q1, Q2, R, testP)
% testP(1) = start point
% testP(2) = step size in point
% testP(3) = number of step

% do GCV .....
GCV=[];
for m=1:testP(3)
    % select time point ....
    V=Av(testP(1)+testP(2)*m,:)';
    k=linspace(-8,-4,50);
    for n=1:length(k)
        lambda=10^k(n);
        [S, L] = sphericalLAP(K,LapK, T, Q1, Q2, R, lambda);
        spV = S*V;
        GCV(n,m)=sum((V-spV).^2)/(1-trace(S)/length(V))^2;
    end
end


[i,ii]= min(GCV);
iii= floor(median(ii));
lambda= 10^k(iii),

% Final S and L .....
[S, L] = sphericalLAP(K,LapK, T, Q1, Q2, R, lambda);

  

