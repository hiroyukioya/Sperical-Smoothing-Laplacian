function [coord] = findSPHE_coord(pid,CH,side)

Aco = xlsread(['~/TT/coordinates/',num2str(pid),'coord.xlsx']);
mni = Aco(:,[1 2 6 7 8]); clear Aco;
    
% Get MNI coordinate..........
MNI=[];
for ni = 1:length(CH)
    f = find(mni(:,1)==CH(ni));
    MNI(ni,:) = mni(f,[3 4 5]);
end
clear ni f 

% read orig.MRI.....
mri = MRIread('~/FS/subjects/icbm152/mri/orig.mgz');

% Apply transform to tkr for AFNI standard mesh
GG=[];
for m=1:size(MNI,1)
    GG(m,:) = [MNI(m,1:3) 0]*inv(mri.vox2ras)'*mri.tkrvox2ras';
end

% Read gifti of AFNI standard surface mesh ......
addpath /home/hiroyuki/gifti-1.8/;
addpath ~/afni_matlab/matlab
if side==1
    g= gifti('~/TT/icbm152/SUMA/std.60.lh.pial.gii');
    gp= gifti('~/TT/icbm152/SUMA/std.60.lh.sphere.gii');
    Bl=importdata('/home/hiroyuki/TT/icbm152/SUMA/std.60.lh.curv.niml.dset');
elseif side==0
    g= gifti('~/TT/icbm152/SUMA/std.60.rh.pial.gii');
    gp= gifti('~/TT/icbm152/SUMA/std.60.rh.sphere.gii');
    Bl=importdata('/home/hiroyuki/TT/icbm152/SUMA/std.60.rh.curv.niml.dset');
end
B = Bl.data;

dist=[]; vind=[];
for k = 1:size(GG,1)
    d = g.vertices - repmat(GG(k,1:3),size(g.vertices,1),1);
    [dist(k),vind(k)] = min(sqrt(sum((d.^2),2)));
end
% Sperical coordinates of contacts .......
coord = gp.vertices(vind,:);
