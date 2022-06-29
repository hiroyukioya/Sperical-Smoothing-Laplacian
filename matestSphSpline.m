function [K, LapK, Q1, Q2, R, T] = matestSphSpline(coord, m)

% Compute matrix required for spherical spline 
% coord = [N channel x (x,y,z)]
% 
%  ref: Carvalhaes and Barros, Int J psychophys (2015)
% _________________________________


% Mapping onto sphere of r=1(m)
az=[]; el=[];
nElec = size(coord,1);
for n=1:size(coord,1)
    [az(n), el(n), r(n)] = cart2sph(coord(n,1),coord(n,2),coord(n,3));
end
ThetaRad = az+pi;
PhiRad = el+pi;
[X,Y,Z] = sph2cart(ThetaRad,PhiRad,1.0);

% Compute cosine distance (Cos(angle)) .......
cos_gamma = zeros(nElec,nElec);    
p = pdist([X' Y' Z'],'euclidean').^2;
pp = squareform(p);
cos_gamma = 1 - pp/2;

% tic
% cos_gamma = zeros(nElec,nElec);    
% % Cosine distance ---------
% % This calculates cos(theta) = 1 - (d^2/2) ......
% for i = 1:nElec
%     for j = 1:nElec      
%       cos_gamma(i,j) = 1 - ( ( (X(i) - X(j))^2 + ...
%         (Y(i) - Y(j))^2 + (Z(i) - Z(j))^2 ) / 2 );
%     end
% end
% toc
r2 = 100;  % to cm ....

G = [];
G0 = 0;
LapG = [];
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
    n = n+1;
end

K = reshape(sum(G,2),nElec,[])/(4*pi);
LapK = reshape(sum(LapG,2),nElec,[])/(4*pi*r2);
T = ones(nElec,1);

% Solve systems of equations.
[Q,R] = qr(T);
R = R(1);
Q1 = Q(:,1);
Q2 = Q(:,2:nElec);




    