function [K, LapK, Q1, Q2, R, T] = matestSphSpline(coord, m)

% Compute matrix required for spherical spline 
% coord = [N channel x (x,y,z)]
%
% _________________________________


% compute radius (r)
r = hypot(coord(:,3), hypot(coord(:,1), coord(:,2)));
r2 = r(1)^2;

% number of electrodes
N = size(coord,1);

% square distance 
sdist = pdist(coord,'euclidean').^                                                                                                              2;
sd2 = squareform(sdist);

% angle 
cos_gamma  = 1 - sd2/(2*r2);

G=[];G0=0;
LapG=[];
epsilon = 1e-10+1;
n=1;
while (1e-10<epsilon)
    Pn = legendre(n, cos_gamma(:));
    a = (2*n+1)/(n*(n+1))^m;
    gm = a*Pn(1,:)';
    G = horzcat(G,gm);
    LapG = horzcat(LapG, -n*(n+1)*gm);
    epsilon = max(abs(G(:,end)-G0));
    G0 = G(:,end);
    n=n+1;
end

K = reshape(sum(G,2),N,[])/(4*pi);
LapK = reshape(sum(LapG,2),N,[])/(4*pi*r2);
T = ones(N,1);

% Solve systems of equations.
[Q,R]= qr(T);
Q1 = Q(:,1);
Q2 = Q(:,2:N);




    