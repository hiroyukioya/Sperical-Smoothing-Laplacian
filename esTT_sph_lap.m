%%  es-TT Spherical Spline Laplacian Potential Analysis ---------

clear; cd ~/TT;
% Load experiment table .......
load ~/TT/HGexp.mat;

% Select runs.......
for  n=1:70
    pid = HGexp(n,1);  
    if HGexp(n,2)<100
            ss=['0',num2str(HGexp(n,2))]
    else
            ss= num2str(HGexp(n,2))
    end  

    [I,A] = load_esTT_pt_info(pid);
    I.explist,  

    analRun=find(strcmp(I.explist,ss)==1);
    eval(['selCH = I.Dpch',num2str(HGexp(n,3)),';']);

    runind = analRun;
    side = HGexp(n,5);

    disp(' ~~~~~~~~~~~~~~~~~ ' )
    disp('  ' )
    fprintf('%s%1.0f%s%s\n','      Patient ',pid,' - Run ',ss);

    if side==0
        disp('           Hemi = RIGHT');
    else
        disp('           Hemi = LEFT');
    end
    disp('  ' )
    disp(' ~~~~~~~~~~~~~~~~~ ' )

    %% Load esTT data .........

    % Read contacts' coordinates ........
    mni=[];
    % [ras, mni,co] = readCoord(pid);
    Aco = xlsread(['~/TT/coordinates/',num2str(pid),'coord.xlsx']);
    mni = Aco(:,[1 2 6 7 8]); clear Aco;

    % read orig.MRI.....
    mri = MRIread('~/FS/subjects/icbm152/mri/orig.mgz');

    % Read surface mesh .......
    [spV, spF] = read_surf('~/FS/subjects/icbm152/surf/lh.sphere');
    [piaV, piaF] = read_surf('~/FS/subjects/icbm152/surf/lh.pial');

    namef = 'sphe'
    basedir = ['~/TT/',num2str(pid),'/',num2str(pid),'-',ss,'-',namef];    

    % Extract NLX data ======
    % Usage :  extract_esTT_data(I,A,analchblk, runid) ............
    I.hpcutoff = 3.0;
    [AvR, Av, CH,  I ] = extract_esTT_data(I,A, selCH, runind);
   
    load(['~/TT/reject',num2str(pid)]);
    eval(['rejch = rej',num2str(pid)]);
    
    % Get MNI coordinate..........
    MNI=[];
    for ni = 1:length(CH)
        f = find(mni(:,1)==CH(ni));
        MNI(ni,:) = mni(f,[3 4 5]);
    end
    clear ni f 
    
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

    % gpreg= gifti('~/TT/icbm152/SUMA/std.60.lh.sphere.reg.gii');

    dist=[]; vind=[];
    for k = 1:size(GG,1)
        d = g.vertices - repmat(GG(k,1:3),size(g.vertices,1),1);
        [dist(k),vind(k)] = min(sqrt(sum((d.^2),2)));
    end
    % Sperical coordinates of contacts .......
    coord = gp.vertices(vind,:);

    addpath ~/es-fMRI
    close all;figure;
    renderFSmeshr_es(gp.vertices, gp.faces+0,B,[-90 0]);%close(figure(1));
    colormap(gray(1026));caxis([-1 2]);hold on;
    % add electrodes.....
    for k=1:size(coord,1);
        plot3(coord(k,1),coord(k,2),coord(k,3),'ks','markersize',22,'markerfacecolor','r');hold on;
    end
    ss= ['~/TT/contacts_on_sphere_',num2str(pid)]
    print(gcf,ss,'-dpng','-r150')
    clear k d m
    ss = ['~/TT/results/data-',num2str(I.pid),'-',I.expid ] ;
    save(ss,'I','coord','g','Av','CH','dist','vind','GG','rejch');
end

%% Spherical Spline Laplacian Analysis ------------------------
clear;
cd ~/TT;
% Load experiment table .......
load ~/TT/HGexp.mat;

for  n=1:70
    pid = HGexp(n,1);  
    if HGexp(n,2)<100
            ss=['0',num2str(HGexp(n,2))]
    else
            ss= num2str(HGexp(n,2))
    end  
    ss = ['~/TT/results/data-',num2str(pid),'-',ss ],
    load(ss)

    close all;
    %  Do spherical spline Laplacian calculations .......
    [K, LapK, Q1, Q2, R, T] = matestSphSpline(coord, 3);

    %  GCV estimate ......
    testP=[2020 2 140];
    [S,L, lambda] = sph_splaplace_FIT(Av, K,LapK, T, Q1, Q2, R, testP);

    % Select time period ....
     f = find(I.tind>-0.5 & I.tind<0.5);
     t = I.tind(f);
     V=Av(f,:);
     % Set bad channels data to zeros ....
     if isempty(rejch)~=1
         fff = ismember(CH,rejch);
         f = find(fff==1)
         V(:,f)=0;
     end
     
     % Laplacian and Spline fit calculation .......
     tLAP=[];spV=[];
     for k= 1:size(V,1)
         tLAP(k,:) = -L*V(k,:)';
         spV(k,:) = S*V(k,:)';
     end
 
    % Compute stats ........
    [zval, thre, sigC] = rms_test(tLAP, 0.01, 0.2, t, I);   % 10-200 ms
    [zval2, thre, sigC] = rms_test(tLAP, 0.1, 0.3, t, I);  % 25-100ms

    % figures ............
%     close all; 
%     f=find(t>0.02); tem = tLAP(f,:);ma = max(abs(tem(:)))
%     figure;
%     FGplot32chN(tLAP(:,1:32), t, I.pid, I.expid ,1 ,[-0.05 0.25], ma*2, CH(1:32));
%     ss=['~/TT/results/sphlap-',num2str(I.pid),'-',I.expid,'-FG']
%     print(gcf,ss,'-dpng','-r150')
% 
%     figure;
%     FGplot32chN(V(:,1:32), t, I.pid, I.expid ,1 ,[-0.05 0.25], 150, CH(1:32))
%     ss=['~/TT/results/CCEP-',num2str(pid),'-',I.expid,'-FG']
%     print(gcf,ss,'-dpng','-r150')
% 
%     figure; 
%     ch=26
%     clf; plot(t,V(:,ch));hold on;plot(t,spV(:,ch),'g'); plot(t,tLAP(:,ch),'r');
%     ylim([-ma*2 ma*2])


    ss = ['~/TT/results/SPHlap',num2str(I.pid),'-',I.expid];
    save(ss, 'I','tLAP','spV','t','CH','lambda','coord','vind','Av','GG','zval','thre','zval2');
end



%%  Response mapping (zval map)
% Load LH surface of afni's standard 60 mesh.
addpath /home/hiroyuki/gifti-1.8/;
addpath ~/afni_matlab/matlab
g= gifti('~/TT/icbm152/SUMA/std.60.lh.pial.gii');
B=importdata('/home/hiroyuki/TT/icbm152/SUMA/std.60.lh.curv.niml.dset');

close all   ;
addpath ~/es-fMRI
[p] = renderFSmeshr_es(g.vertices, g.faces+0,B.data,[-90 0]);
o=findobj(gca,'type','light');colormap(gray(1026))
delete(o)
light('position',[-6 6 10],'color','w')
light('position',[-6 2 10],'color', 'w')
light('position',[0 12 -10],'color','w')
light('position',[-5 -5 -5],'color','w')
light('position',[-20 -22 10],'color','w')
light('position',[10 0 10],'color','w')
light('position',[1 20 0],'color','w')
p=findobj(gcf,'type','patch');   
set(p,'edgealpha',0,'specularexponent',0.9,'specularstrength',0.25,...
    'specularColorReflectance',0.9,'AmbientStrength',1,'Diffusestrength',0.8);
set(gcf,'color','w'); 
set(gca,'xcolor','k','ycolor','k','fontsize',24)
close(figure(1)); hold on;
 
% ss=['~/TT/results/MNI_Lh_Lateral_Surface'];
% set(gca,'fontsize',24); ax=axis;
% print(gcf,ss,'-dpng','-r300');
% saveas(gcf,ss)


%% //////   HK smoothing of Z-vals   ////////

addpath ~/Tone/HGmap
cd ~/TT;
% Load experiment table ................................
load ~/TT/HGexp.mat;
exptoanal = [1:70];
GGG=[]; Zv = [];
sdata = zeros(size(g.vertices,1),length(exptoanal));sdata2=sdata;
for n = 1:length(exptoanal)
    n,
    pid = HGexp(n,1);
    if HGexp(n,2)<100
            s=['0',num2str(HGexp(n,2))];
    else
            s= num2str(HGexp(n,2));
    end  
    ss = ['~/TT/results/SPHlap',num2str(pid),'-',s];
    disp(' ~~~~~~~~~~~~~~~~~ ' )
    disp('  ' )
    fprintf('%s%s\n','   Loading .....  ',ss);
    load(ss);
    
    f = find(GG(:,1)>=0);
    GG(f,1) = -GG(f,1);
    % Find nearest vertex .................................
    dist = []; vind = [];
    for k = 1:size(GG,1)
        d = g.vertices - repmat(GG(k,1:3),size(g.vertices,1),1);
        [dist(k),vind(k)] = min(sqrt(sum((d.^2),2)));
    end
%     plot(dist);waitforbuttonpress;
    Z = zeros(size(g.vertices,1),1);
%     Z2 = zeros(size(g.vertices,1),1);
    Z(vind) = zval;
%     Z2(vind) = zval2;
%     [sdata(:,n)] = heatKernelSmoothZ(g.vertices, g.faces, Z, 1.75, 1);
    [out] = hk_smooth(Z, g, 3.3,1);   % Smoothing ... (param 3.3 = 6mm FWHM)
%     [out2] = hk_smooth(Z2, g, 3.3,1);   % Smoothing ... (param 3.3 = 6mm FWHM)
    
    R = max(Z)/max(out);
    sdata(:,n) = out*R;
    
%     R = max(Z2)/max(out2);
%     sdata2(:,n) = out*R;
end
ss=['~/TT/results/sphe_spline_lap_HKsdata'];
save(ss,'sdata','g');


% Surface anatomy image ----
close all
renderFSmeshr_es(g.vertices, g.faces+0, B.data ,[-90 0]);
close(figure(1)); colormap(gray);
ss=['~/TT/results/MNI152_Frontal_Lateral'];
set(gca,'fontsize',24);
set(gca,'xcolor','w','ycolor','w','fontsize',24); grid off;

print(gcf,ss,'-dpng','-r300');
saveas(gcf,ss)
 
% zval mapping ===========================
close all;
ff=[];
f = find(HGexp(:,4)==1);
ff = find(HGexp(:,1)==405 | HGexp(:,1)==394); % Bad subject removal
f = setdiff(f,ff);

ms = mean(sdata(:,f),2);
renderFSmeshr_es(g.vertices, g.faces+0,ms  ,[-90 0]);
cmap=flipud(colormap(hot(128))); 
cmap = cmap([1:75],:);
colormap(cmap); caxis([3 10]); o=findobj(gca,'type','light');delete(o)
light('position',[-6 6 0],'color','w');
close(figure(1));
set(gca,'xcolor','w','ycolor','w','fontsize',24); grid off;


ss=['~/TT/results/MNI_NONCore_Zval_Frontal_Lateral'];
set(gca,'fontsize',24)
print(gcf,ss,'-dpng','-r300');
saveas(gcf,ss)
 




% kernel weight plot............
K=@(x,sigma) exp(-x./(4*sigma))./sum(exp(-x./(4*sigma)));
figure;
x=K([0:.1:200],3.3);  ma = max(x); plot(sqrt([0:.1:200]), x/ma,'linewidth',2);
ylim([0 1]); grid on;xlabel('distance (mm)'); ylabel('normalized weight')
xlim([0 10])
ss=['~/TT/heatkernel-weight-plot'];
print(ss,'-dpng','-r150');
 
 
%%   Contacts mapping .........
% clear; 
cd ~/TT;
load ~/TT/HGexp.mat;
mri = MRIread('~/FS/subjects/icbm152/mri/orig.mgz');
ptlist = [369;372;376; 384;394;399;405;413;423;429;460];
L = size(HGexp,1);
rz = zeros(32,4,L);
stimcoord=[];
for n=1:L
    pid = HGexp(n,1);
    if HGexp(n,2)<100
            ss=['~/TT/results/SPHlap',num2str(HGexp(n,1)),'-0',num2str(HGexp(n,2))]
    else
            ss=['~/TT/results/SPHlap',num2str(HGexp(n,1)),'-',num2str(HGexp(n,2))]
    end  
    
    % Load data ......
    load(ss);
    
    mni = [];
    Aco = xlsread(['~/TT/coordinates/',num2str(pid),'coord.xlsx']);
    mni = Aco(:,[1 2 6 7 8]); clear Aco;
    
    % Apply transform to tkr for AFNI standard mesh
    MNI=[];
    for m=1:size(mni,1)
        MNI(m,:) = [mni(m,3:5) 0]*inv(mri.vox2ras)'*mri.tkrvox2ras';
    end
    
    % Apply transform to tkr for FS
    MNIfs=[];
    for m=1:size(mni,1)
        MNIfs(m,:) = [mni(m,3:5) 1]*inv(mri.vox2ras)'*mri.tkrvox2ras';
    end

    % Find Stimulated channel coordinates .....
    iii = find(ismember(mni(:,1),I.stimch));
    stimcoord(n,:) = (MNI(iii(1),[1 2 3]) + MNI(iii(2),[1 2 3]))/2; % For afni std mesh
    
    % Find Stimulated channel coordinates .....
    iii = find(ismember(mni(:,1),I.stimch));
    stimcoordFS(n,:) = (MNIfs(iii(1),[1 2 3]) + MNIfs(iii(2),[1 2 3]))/2; % For afni std mesh

    % Find Recorded channel coordinates .....
    for kk  = 1:length(CH)
          ff = find(mni(:,1)==CH(kk));
          rz(kk,:,n) = [MNI(ff,[1 2 3]) zval(kk)'];
    end
    
    % Find Recorded channel coordinates (for FS mesh)
    for kk  = 1:length(CH)
          ff = find(mni(:,1)==CH(kk));
          rzfs(kk,:,n) = [MNIfs(ff,[1 2 3]) zval(kk)'];
    end
end

clear n L m f MNI iii Av AvR CH co I inp ras mni ss tlap thre ans f 
%%  Find nearest vertex ......................
clear zv
% [IMPORTANT] Use "standard mesh" created by AFNI ..........
A=csvread('/home/hiroyuki/TT/std_icbm152_lh60_pial.csv',0,0);
ind = A(1,1);
std60Vert = A([2:ind+1],[1 2 3]);  size(std60Vert)
std60Face = A([ind+2:end],[1 2 3]); size(std60Face)
B=csvread('/home/hiroyuki/TT/std.60.lh.curv.csv',0,0);

j = 0; RS=[];
for n = 1:size(HGexp,1)
    for k = 1:32
        j = j+1;
        
        coord= rz(k,[1:3],n);
        if coord(1)>= 0
            recside= 0;
        else 
            recside = 1;
        end
        co = rz(k,[1 2 3],n);
        if recside==0
            co(1) = -co(1);
        end
        d = std60Vert - co;
        [dist(j),ii] = min(sqrt(sum((d.^2),2)));
        
        RS(k,4,n) = ii;
        RS(k,3,n) = recside;
        RS(k,1,n) = HGexp(n,1);
        RS(k,2,n) = HGexp(n,2);
        RS(k,[5 6 7],n) = rz(k,[1 2 3],n); % For AFNI's std mesh
        tem =stimcoordFS(n,[1 2 3]);
        if abs(tem(1))<33
            tem(1) =35*sign(tem(1));
        end
        RS(k,[8 9 10],n) = tem;
        RS(k,[11 12 13],n) = rzfs(k,[1 2 3],n); % For FS's mesh
        RS(k,14,n) = rz(k,4,n);
    end
end
% Combined data = RS [contacts x data x runs]
% RS(:,2,:) =>
%   [pid - run# - recSide - vertexid - RecordMNI-AFNI(5-7) - ...
%             StimMNI-FS(8-10) -  RecordMNI-FS(11 12 13)-  Zval]       

clear n k j d ii recside stimmni zval pt ;

%%

%///////////   Contacts mapping  /////////////%
addpath /home/hiroyuki/gifti-1.8/;
addpath ~/afni_matlab/matlab;
addpath   /home/hiroyuki/es-fMRI
g= gifti('~/TT/icbm152/SUMA/std.60.lh.pial.gii');
B=importdata('/home/hiroyuki/TT/icbm152/SUMA/std.60.lh.curv.niml.dset');

close all;renderFSmeshr_es(g.vertices, g.faces+0,B.data,[-90 0]);
o=findobj(gca,'type','light');colormap(gray(1026))
set(gcf,'color','w'); 
set(gca,'xcolor','w','ycolor','w','fontsize',24)
close(figure(1)); hold on;

ss=['~/TT/results/MNI_Lh_Lateral_Surface'];
set(gca,'fontsize',24); ax=axis;
print(gcf,ss,'-dpng','-r300');
saveas(gcf,ss)

i=[]; ii=[];
o=findobj(gca,'type','line');delete(o)
 
J = find(HGexp(:,4)<=3);
M=squeeze(RS(:, end, J)); 
% d= [1:8 9 16 17 24 25:32];   M(d,:,:)=0;
f=find(M>=0 );
[i,ii]=ind2sub([size(M)],f);
FGcoord=[];

pt = unique(HGexp(:,1));
for n=1:length(pt)
     f = find(HGexp(:,1)==pt(n));
     f=f(1);
     tem=RS(:, [5 6 7], f);
     ff=find(tem(:,1)>0);
     tem(ff,1) = -tem(ff,1);
     FGcoord(:,:,n)=tem;
end

pp = [1 2 3 4 6 8 9 10 11]
for n=1:length(pp)
    nn = pp(n);
    hold on;plot3(FGcoord(:,1,nn)-10,FGcoord(:,2,nn),FGcoord(:,3,nn),'bo','linewidth',4,'markersize',20);
end

ss=['~/TT/results/Contacts_mapping_Frontal_grids'];
set(gca,'fontsize',24)
print(gcf,ss,'-dpng','-r300');


%%   STP sites ...............................
addpath /home/hiroyuki/MRI_tone/HGmap
addpath ~/Tone/HGmap:~/es-fMRI  
addpath  /usr/local/FS6.0/matlab;

% Load subject's STP patch ......  
close all;
[mypatch] = load_PT_patch_es('icbm152','lh');

% Draw STP surface .....
renderFSmeshr_es(mypatch.vertsPia, mypatch.faces, mypatch.curv ,...
    [180 65]);

p=findobj(gcf,'type','patch');   
set(p,'edgealpha',0,'specularexponent',0.9,'specularstrength',0.25,...
    'specularColorReflectance',0.9,'AmbientStrength',1,'Diffusestrength',0.8);
set(gcf,'color','w'); 
set(gca,'xcolor','k','ycolor','k','fontsize',24)
colormap(gray);
o=findobj(gca,'type','light');colormap(gray(1026))
delete(o)
light('position',[50 6 100],'color',[0.6 0.6 0.6]);
light('position',[-60 2 100],'color', [0.6 0.6 0.6]);
light('position',[50 60 100],'color',[0.6 0.6 0.6]);
light('position',[-60 -60 100],'color', [0.6 0.6 0.6]);
light('position',[0 40 10],'color', 'w');
ss=['~/TT/results/MNI_lhSTP_surface'];
set(gca,'fontsize',24);set(gca,'fontsize',24);set(gca,'xcolor','w','ycolor','w','fontsize',24)

set(gcf,'paperunit','inches','papersize',[3 3],'InvertHardcopy','off');
print(gcf,ss,'-dpng','-r300');
saveas(gcf,ss)
ax=axis;

% Stimulating contacts ......
stC = squeeze(RS(1,[8 9 10],:))';
roi = HGexp(:,4);

% Plot contacts .....
o=findobj(gca,'type','line');delete(o)
o=findobj(gca,'type','text');delete(o)

for n=1:size(HGexp,1)
     tem = stC(n,:);
    if mean(tem(1))>0
        tem(1) = -tem(1);
    end
    if roi(n)==1
        plot3(tem(1),tem(2),tem(3)+0,'o','markeredgecolor','b','markersize',20,...
        'linewidth',4);
     text(tem(1),tem(2),tem(3)+1,num2str(n),'fontsize',24)
    elseif roi(n)==2
         plot3(tem(1),tem(2),tem(3)+0,'o','markeredgecolor','r','markersize',20,...
          'linewidth',4);
      text(tem(1),tem(2),tem(3)+1,num2str(n),'fontsize',24)
    elseif roi(n)==3
         plot3(tem(1),tem(2),tem(3)+0,'o','markeredgecolor','r','markersize',20,...
          'linewidth',4);
       text(tem(1),tem(2),tem(3)+1,num2str(n),'fontsize',24)
    end
end

ss=['~/TT/results/MNI_lhSTP_surface_w_contacts'];
set(gca,'fontsize',24);set(gca,'xcolor','w','ycolor','w','fontsize',24)
ax=axis


% set(gcf,'paperunit','inches','papersize',[3 3],'InvertHardcopy','off');axis(ax)
o=findobj(gca,'type','patch');delete(o); axis(ax);

%%
%%  ROI average waveform comparisons......
clear; cd ~/TT;
load ~/TT/HGexp.mat;
ptlist = [369;372;376; 384;394;399;405;413;423;429;460];
L = size(HGexp,1);
MAV = zeros(1000,L);MAV2=MAV;
Z=[];
for n = 1:L
    if HGexp(n,2)<100
            ss = ['~/TT/results/SPHlap',num2str(HGexp(n,1)),'-0',num2str(HGexp(n,2))];
    else
            ss = ['~/TT/results/SPHlap',num2str(HGexp(n,1)),'-',num2str(HGexp(n,2))];           
    end 
    % Load data ......
    load(ss);
    
    % High-pass filtering ...............
%     tLAPhp=[];
%     for nc=1:size(tLAP,2)
%         [tLAPhp(:,nc)] = IIRHPfilter(tLAP(:,nc), I.fs, 15);
%     end
    tLAPhp = tLAP;
    
    %Normalize tlap waveform................
    f = find(t>=0.01 & t<= 0.5);
    tem = sqrt(mean(tLAPhp(f,:).^2,1));
%      tem = max((tLAP(f,:)),[],1);
    nlap = tLAPhp./repmat(tem,size(tLAPhp,1),1);
    
    if HGexp(n,1)==369
        [i,ifgT] = intersect(CH,[227:230 236:238]);
    elseif HGexp(n,1)==372
        [i,ifgT] = intersect(CH,[241:246 249:253]);
    elseif HGexp(n,1)==376
        [i,ifgT] = intersect(CH,[130:133 141]);
    elseif HGexp(n,1)==384
        [i,ifgT] = intersect(CH,[]);
    elseif HGexp(n,1)==394
        [i,ifgT] = intersect(CH,[65:68 73:76]);
         ifgT=[];  % No IFG contact
    elseif HGexp(n,1)==399
        [i,ifgT] = intersect(CH,[131:134 140:142]);
    elseif HGexp(n,1)==405
        [i,ifgT] = intersect(CH,[]);
         ifgT=[];  % No IFG contact
    elseif HGexp(n,1)==413
        [i,ifgT] = intersect(CH,[148:151 158:159]);
    elseif HGexp(n,1)==423
        [i,ifgT] = intersect(CH,[187:190 180:183 175]);
    elseif HGexp(n,1)==429
        [i,ifgT] = intersect(CH,[179:182 190]);
    elseif HGexp(n,1)==460
        [i,ifgT] = intersect(CH,[219:222 211:213]);
    end
    
    ii = find(zval(ifgT)>3);
    ifg = ifgT(ii);
    clf; 
    if isempty(ifg)~=1
        MAV(1:size(tLAP,1),n) = mean(nlap(:,ifg), 2);
    else 
        MAV(1:size(tLAP,1),n) = zeros(2000,1);
    end
    
    [k,kk] = max(zval(ifgT));
    if isempty(kk)~=1
        MAV2(1:size(tLAP,1),n) = nlap(:,ifgT(kk));    
    else
        MAV2(1:size(tLAP,1),n) = 0;    
    end
end

f=find(t>0.01);
ff= find(mean(MAV,1)==0);

core = setdiff(find(HGexp(:,4)==1),ff);
noncore = setdiff(find(HGexp(:,4)>=2),ff);

% BASELINE   correction.
f = find(t>-0.05 & t<=-0.01);
m = mean(MAV(f,:),1);
MAV = MAV-repmat(m,size(MAV,1),1);
% delete bad waveforms
core([ 17 18 20 22])=[];

% MAV(:,core([1 5 8 11 12 16]))=-MAV(:,core([1 5 8 11 12 16]));

K1 = mean(-MAV(:,noncore),2);
K2 = mean(-MAV(:,core),2);

% draw waveforms with se 
addpath /home/hiroyuki/MRI_tone;
% close all;
figure;
plotRPwSE(K1, std(MAV(:,noncore),[],2)/sqrt(length(noncore)),...
    1, t,'r');hold on;
plotRPwSE(K2, std(MAV(:,core),[],2)/sqrt(length(core)), ...
    1, t,'b')
xlim([-0.05 0.25]);ylim([-1.5 1.5]);
grid on; xlabel('Time (s)'); ylabel('Normalized laplacian amplituude (\muV/cm^2)')
set(gca,'gridlinestyle','-','fontsize',26)
ss=['~/TT/HG-IFG-overall-waveforms'];
set(gca,'linewidth',1,'tickdir','out');
legend('','Lateral HG','','Medial HG')

print(gcf,ss,'-dpng','-r300');
saveas(gcf,ss)


% clf;
plot(t, mean(MAV(:,core),2),'linewidth',1.5,'color','b');
xlim([-0.05 0.25]);ylim([-ma*1.5 ma*1.5]);hold on;
plot(t, mean(MAV(:,noncore),2),'linewidth',1.5,'color','r');
xlim([-0.05 0.25]);ylim([-ma*1.5 ma*1.5]);
grid on; xlabel('Time (s)'); ylabel('Normalized  \muV/cm^2')
set(gca,'gridlinestyle','-','fontsize',24)
ss=['~/TT/HG-IFG-overall-waveform-comparison'];
print(gcf,ss,'-dpng','-r300');

%%  HIPPO ROI average waveform comparisons......
clear; cd ~/TT;
load('~/TT/Hippoexp.mat');Hippoexp= [Hippoexp Hippoexp(:,3) Hippoexp(:,3)];
HGexp=Hippoexp;
selrun = [1];

MAV = zeros(1000,L);MAV2=MAV;
Z=[];
for n = 1:length(selrun)
    m = selrun(n);
    if HGexp(n,2)<100
            ss = ['~/TT/results/SPHlap',num2str(HGexp(m,1)),'-0',num2str(HGexp(m,2))];
    else
            ss = ['~/TT/results/SPHlap',num2str(HGexp(m,1)),'-',num2str(HGexp(m,2))];           
    end 
    % Load data ......
    load(ss);
    
    % High-pass filtering ...............
%     tLAPhp=[];
%     for nc=1:size(tLAP,2)
%         [tLAPhp(:,nc)] = IIRHPfilter(tLAP(:,nc), I.fs, 15);
%     end
    tLAPhp = tLAP;
    
    %Normalize tlap waveform................
    f = find(t>=0.01 & t<= 0.5);
    tem = sqrt(mean(tLAPhp(f,:).^2,1));
%      tem = max((tLAP(f,:)),[],1);
    nlap = tLAPhp./repmat(tem,size(tLAPhp,1),1);
    
    if HGexp(n,1)==369
        [i,ifgT] = intersect(CH,[227:230 236:238]);
    elseif HGexp(n,1)==372
        [i,ifgT] = intersect(CH,[241:246 249:253]);
    elseif HGexp(n,1)==376
        [i,ifgT] = intersect(CH,[130:133 141]);
    elseif HGexp(n,1)==384
        [i,ifgT] = intersect(CH,[]);
    elseif HGexp(n,1)==394
        [i,ifgT] = intersect(CH,[65:68 73:76]);
         ifgT=[];  % No IFG contact
    elseif HGexp(n,1)==399
        [i,ifgT] = intersect(CH,[131:134 140:142]);
    elseif HGexp(n,1)==405
        [i,ifgT] = intersect(CH,[]);
         ifgT=[];  % No IFG contact
    elseif HGexp(n,1)==413
        [i,ifgT] = intersect(CH,[148:151 158:159]);
    elseif HGexp(n,1)==423
        [i,ifgT] = intersect(CH,[187:190 180:183 175]);
    elseif HGexp(n,1)==429
        [i,ifgT] = intersect(CH,[179:182 190]);
    elseif HGexp(n,1)==460
        [i,ifgT] = intersect(CH,[219:222 211:213]);
    end
    
    ii = find(zval(ifgT)>3);
    ifg = ifgT(ii);
    clf; 
    if isempty(ifg)~=1
        MAV(1:size(tLAP,1),n) = mean(nlap(:,ifg), 2);
    else 
        MAV(1:size(tLAP,1),n) = zeros(2000,1);
    end
    
    [k,kk] = max(zval(ifgT));
    if isempty(kk)~=1
        MAV2(1:size(tLAP,1),n) = nlap(:,ifgT(kk));    
    else
        MAV2(1:size(tLAP,1),n) = 0;    
    end
end

f=find(t>0.01);
ff= find(mean(MAV,1)==0);

core = setdiff(find(HGexp(:,4)==1),ff);
noncore = setdiff(find(HGexp(:,4)>=2),ff);

% BASELINE   correction.
f = find(t>-0.05 & t<=-0.01);
m = mean(MAV(f,:),1);
MAV = MAV-repmat(m,size(MAV,1),1);
% delete bad waveforms
core([ 17 18 20 22])=[];

% MAV(:,core([1 5 8 11 12 16]))=-MAV(:,core([1 5 8 11 12 16]));

K1 = mean(-MAV(:,noncore),2);
K2 = mean(-MAV(:,core),2);

% draw waveforms with se 
addpath /home/hiroyuki/MRI_tone;
% close all;
figure;
plotRPwSE(K1, std(MAV(:,noncore),[],2)/sqrt(length(noncore)),...
    1, t,'r');hold on;
plotRPwSE(K2, std(MAV(:,core),[],2)/sqrt(length(core)), ...
    1, t,'b')
xlim([-0.05 0.25]);ylim([-1.5 1.5]);
grid on; xlabel('Time (s)'); ylabel('Normalized laplacian amplituude (\muV/cm^2)')
set(gca,'gridlinestyle','-','fontsize',26)
ss=['~/TT/HG-IFG-overall-waveforms'];
set(gca,'linewidth',1,'tickdir','out');
legend('','Lateral HG','','Medial HG')

print(gcf,ss,'-dpng','-r300');
saveas(gcf,ss)



%%  Latency measurements ....
fid = fopen('~/TT/latency.csv','w');
for n=1:70
    plot(t,MAV2(:,n)); xlim([-0.05 0.15]);
    [x,y] = ginput(1);
    x=x*1000;
%     fprintf(fid,'%3.1f\t%3.1f\t%3.1f\t%3.1f\n',x(1),x(2),x(3),x(4));
    fprintf(fid,'%3.1f\n',x(1));
end
fclose(fid);
 

fl = csvread('~/TT/latency.csv')
f= find(fl>0);
core = find(HGexp(:,4)==1);
noncore = find(HGexp(:,4)>=2);
c1 = intersect(core,f);
c2 = intersect(noncore,f);

mean(fl(c1)), std(fl(c1))
mean(fl(c2)), std(fl(c2))

[H,P]=ttest2(fl(c1),fl(c2))

close all;
 for n=1:length(core)
    plot(t,-MAV(:,core(n)));xlim([-0.05 0.2]);ylim([-5 5]); title(num2str(n))
    ginput(1);
end

%% |||||||      High-gamma power movie creation     |||||||||
clear n sel* mni Aco pia* 
% HG stimulation esTT....
load('~/TT/HGexp.mat');

% if Hippo stimulation es-TT....
load('~/TT/Hippoexp.mat');Hippoexp= [Hippoexp Hippoexp(:,3) Hippoexp(:,3)];
HGexp=Hippoexp;

% mode 1= high gamma,    mode 2 = Voltage Measurements.
mode = 2;

for  n = 9:10
    pid = HGexp(n,1);  
    if HGexp(n,2)<100
            ss=['0',num2str(HGexp(n,2))];
    else
            ss= num2str(HGexp(n,2));
    end  

    [I,A] = load_esTT_pt_info(pid);
    I.explist;

    analRun=find(strcmp(I.explist,ss)==1);
    load(['~/TT/LateralGrid',num2str(pid),'.mat'])

    runind = analRun;
    side = HGexp(n,5);

    disp(' ~~~~~~~~~~~~~~~~~ ' )
    disp('  ' )
    fprintf('%1.0f%s%1.0f%s%s\n', n,'      Patient ',pid,' - Run ',ss);
    
    % Read contacts' coordinates ........
    mni=[];
    % [ras, mni,co] = readCoord(pid);
    Aco = xlsread(['~/TT/coordinates/',num2str(pid),'coord.xlsx']);
    mni = Aco(:,[1 2 6 7 8]); clear Aco;

    % read orig.MRI.....
    mri = MRIread('~/FS/subjects/icbm152/mri/orig.mgz');

    % Read surface mesh .......
    [piaV, piaF] = read_surf('~/FS/subjects/icbm152/surf/lh.pial');

    % Extract NLX data =============
    I.hpcutoff = 6.0;  % set HP cutoff freq. 
    [AvR, Av, CH,  I , ndat, gsig, hsig] = extract_esTT_data(I,A, selCH, runind);
        
    % Band pass filtering ============
    f = find(I.tind>=-0.6 & I.tind<=0.6);
    t = I.tind(f);
    baseline = find(t>=-0.5 & t<=-0.2);
    erbp = zeros(length(t),size(ndat,2));
    
    % Calc ERBP =================+
    parfor c = 1: size(ndat,3)
        tem = abs(hilbert(WindowedBandPassFilter(ndat(f,:,c),I.fs,128,70,140)));
        b = median(tem(baseline,:),1);
        seltrial = cell2mat(I.tr(c));
        b=b(seltrial);  b = repmat(b,size(tem,1),1); 
        erbp(:,c) = mean(20*log10(tem(:,seltrial)./b),2);
    end
        bia = mean(mean(erbp(200:800,:),1));
        erbp =erbp-bia;
        
    % Get MNI coordinate  ............
    MNI=[];
    for ni = 1:length(CH)
        f = find(mni(:,1)==CH(ni));
        MNI(ni,:) = mni(f,[3 4 5]);
    end
    clear ni f bia b seltrial tem
    
    % Apply transform to tkr for AFNI standard mesh ..............
    GG=[];
    for m=1:size(MNI,1)
        GG(m,:) = [MNI(m,1:3) 0]*inv(mri.vox2ras)'*mri.tkrvox2ras';
    end
    
    % load mesh ........
    addpath /home/hiroyuki/gifti-1.8/;
    addpath ~/afni_matlab/matlab
    if side==1
        g= gifti('~/TT/icbm152/SUMA/std.60.lh.pial.gii');
        Bl=importdata('/home/hiroyuki/TT/icbm152/SUMA/std.60.lh.curv.niml.dset');
    elseif side==0
        g= gifti('~/TT/icbm152/SUMA/std.60.rh.pial.gii');
        gp= gifti('~/TT/icbm152/SUMA/std.60.rh.sphere.gii');
        Bl=importdata('/home/hiroyuki/TT/icbm152/SUMA/std.60.rh.curv.niml.dset');
    end
    B = Bl.data;
    
    %////////////////////////////////////////////////////////////%
    % Draw Brain Surface 
    addpath ~/es-fMRI
    close all;figure;
    renderFSmeshr_es(g.vertices, g.faces+0,B,[-90 0]);close(figure(1));
    colormap(gray(1026));caxis([-0.8 2]);H=findobj(gcf,'type','patch');   
    set(H,'edgealpha',0,'specularexponent',0.35,'specularstrength',0.2,...
    'specularColorReflectance',1,'AmbientStrength',1,'Diffusestrength',0.5);
    set(gcf,'color','w'); set(gca,'xcolor','w','ycolor','w','fontsize',24)
    if side ==0 
        view(90,0);
        o=findobj(gca,'type','light');colormap(gray(1024));delete(o)
        light('position',[1000 0 0],'color','w');
        light('position',[100 150 0],'color','w');
        light('position',[100 -150 0],'color','w');
    else
        view(-90,0);
        o=findobj(gca,'type','light');colormap(gray(1024));delete(o)
        light('position',[-1000 0 0],'color','w');
        light('position',[-100 150 0],'color','w');
         light('position',[-100 -150 0],'color','w');
    end
    %  //////////////////////////////////////////////////////// %
       
    % Time bins ......
    tperiod = 0.000:0.002:0.4;
    colscale = parula(256);  % define colorscale....
    a=colorcet('CBD2');
    colscale=a;
    maxz = 2;   % max/min colorscale value....
    vv = floor(length(colscale)/(2*maxz)); 
    
    % < Create AVI video file >
    % Use ERBP ------------------------------------------
    if mode ==1
        ss = ['~/TT/erbp/movie_estt_',num2str(pid),'_',I.expid,'.avi'];
        v = VideoWriter(ss,'Motion JPEG AVI');
        v.FrameRate = 8; 
        v.Quality = 50;
        open(v)
        for ii = 1:length(tperiod)-1
            dur = find(t>=tperiod(ii) & t<tperiod(ii+1));
            ctime = floor(1000*(tperiod(ii)+tperiod(ii+1))/2);
            
            % Thresholding .......
            f00=find(t>-0.5 & t<-0.2);
            x00=erbp(f00,:);x00=x00(:);
            [lim]=prctile(x00,[10 90]);
            f000=find(erbp>lim(1) & erbp<lim(2));
            erbp(f000)=0;

            if ctime < 0.008*1000
                x = zeros(length(CH),1);
            else
                x = mean(erbp(dur,:),1);
            end
            title([num2str(ctime) ' ms'],'fontsize',36);
            o=findobj(gca,'type','line');delete(o);

            for k=1:size(GG,1)
                if x(k)>=0
                    cc = floor(min(x(k)*vv,128))+128;
                else
                    cc = floor(max(x(k)*vv,-128))+129;
                end
                if side==0
                    h2=plot3(GG(k,1)+20, GG(k,2),GG(k,3),'o','markersize',28,...
                        'markerfacecolor',[colscale(cc,:)],'markeredgecolor',[colscale(cc,:)]); 
                    hold on; 
                else 
                    h2=plot3(GG(k,1)-20, GG(k,2),GG(k,3),'o','markersize',28,...
                        'markerfacecolor',[colscale(cc,:)],'markeredgecolor',[colscale(cc,:)]); 
                    hold on; 
                end
            end      

            % Set colorbar property .....
            if ii==1
                yy=[];
                o=get(gcf,'ColorMap'); nl = size(o,1);
                yy(1:nl,:)=o; yy(nl+1:nl+256,:)=colscale;
                set(gcf,'colormap',yy); ca = caxis;H=linspace(ca(1), ca(2), size(yy,1));
                h=colorbar;  
                set(h,'Location','south','position',[0.4 0.15 0.22 0.02],'Limits',[H(1025) H(end)],...
                    'ticks',[H(1025) H(end)],'ticklabels',[-maxz maxz],'fontsize',36);
            end
            drawnow;
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        close(v);
        clear k cc x o dur ctime ans H h f ff i ii ni cc o h2 yy o H 
        ss = ['~/TT/erbp/erbp_estt_',num2str(pid),'_',I.expid,'.mat']

    % Use Voltage waveform ------------------------
    elseif mode==2
        clear m n tr ttr ;
        % define colormap (should be 256 bins)
        a=colorcet('CBD2');
        colscale=a;
%         colscale = parula(256);
%         colscale(128-10:128+10,:)=repmat([0.7 0.7 0.7],21,1);
        gsm = []; erbp = [];
        f = find(I.tind>=-0.6 & I.tind<=0.6);
        t = I.tind(f);
        baseline = find(t>=-0.5 & t<=-0.2);
        
        % High-pass filtering ...............
        ndath = ndat;
        
        % Common-average re-reference ...........
        gs= zeros(size(ndat));
        for tr = 1:size(ndat,2)
            m = mean(ndath(:,tr,:),3);
            gs(:,tr,:) = ndath(:,tr,:) - repmat(m,1,1,size(ndath,3));
        end
        
        % Mean gs with trial rejection .........
        for u = 1: size(ndath,3)
            ttr = cell2mat(I.tr(u));
            gsm(:,u) = mean(gs(:,ttr,u),2);
        end
        clear u ttr ;  
        
        % GS= Induced potentials (applied mean subtraction) .........
        GS=[];
        GS = gs - permute(repmat(gsm,1,1,size(gs,2)),[1 3 2]);
        
        % Mean abs((GS)) with trial rejection ..........
        for u = 1: size(ndath,3)
            ttr = cell2mat(I.tr(u));
            mGS(:,u) = mean(abs(GS(:,ttr,u)),2);
        end
         
        % Normalization and calculate dB change  .........
        erbp = mGS(f,:);
        baseMe_erbp = repmat(median(erbp(baseline,:),1),size(erbp,1),1);        
%         baseSD_erbp = repmat(std(erbp(baseline,:),1),size(erbp,1),1);
        erbp = 20*log10(erbp./baseMe_erbp);  % dB scale
        
        % Apply threshold ........
        r = -0.55 + (-0.075+0.5)*rand(1000,1);
        for k = 1:1000
            d= find(t>=r(k) & t<r(k)+0.005);
            basem(k,:) =mean(erbp(d,:),1);
        end
        TH= prctile(basem,[2.5 97.5],1);
        
        erbpth = erbp;
        for k=1:size(erbp,2)
            i = find(erbp(:,k)>=TH(1,k) & erbp(:,k)<=TH(2,k));
            erbpth(i,k)=0;
        end
         
        ss = ['~/TT/erbp/AEP_movie_estt_',num2str(pid),'_',I.expid,'.avi'];
        v = VideoWriter(ss,'Motion JPEG AVI');
        v.FrameRate = 8; 
        v.Quality = 50;
        maxz = 3;   % max/min colorscale value....
        tperiod = -0.0075:0.005:0.5;
        
        open(v)
        for ii = 1:length(tperiod)-1
            dur = find(t>=tperiod(ii) & t<tperiod(ii+1));
            ctime = floor(1000*(tperiod(ii)+tperiod(ii+1))/2);         

            if ctime < 0.004*1000
                x = zeros(length(CH),1);
            else
                x = mean(erbpth(dur,:),1);
            end
            
            title([num2str(ctime) ' ms'],'fontsize',36);
            o=findobj(gca,'type','line');delete(o);

            for k=1:size(GG,1)
                if x(k)>=0
                    cc = floor(min(x(k)*vv,128))+128;
                else
                    cc = floor(max(x(k)*vv,-128))+129;
                end
                if side==0
                    h2=plot3(GG(k,1)+20, GG(k,2),GG(k,3),'o','markersize',28,...
                        'markerfacecolor',[colscale(cc,:)],'markeredgecolor',[colscale(cc,:)]); 
                    hold on; 
                else 
                    h2=plot3(GG(k,1)-20, GG(k,2),GG(k,3),'o','markersize',28,...
                        'markerfacecolor',[colscale(cc,:)],'markeredgecolor',[colscale(cc,:)]); 
                    hold on; 
                end
            end      

            % Set colorbar property .....
            if ii==1
                yy=[];
                o=get(gcf,'ColorMap'); nl = size(o,1);
                yy(1:nl,:)=o; yy(nl+1:nl+256,:)=colscale;
                set(gcf,'colormap',yy); ca = caxis;H=linspace(ca(1), ca(2), size(yy,1));
                h=colorbar;  
                set(h,'Location','south','position',[0.4 0.15 0.22 0.02],'Limits',[H(1025) H(end)],...
                    'ticks',[H(1025) H(end)],'ticklabels',[-maxz maxz],'fontsize',36);
            end
            drawnow;
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        close(v);
        clear k cc x o dur ctime ans H h f ff i ii ni cc o h2 yy o H 
        ss = ['~/TT/erbp/induce_estt_',num2str(pid),'_',I.expid,'.mat']
    end
    save(ss,'erbp','erbpth','GG','side','t','CH','Av','I','mni','TH')

end

%% -------------   Make averaged movie  --------------
clear; cd ~/TT;
hipp = 1;

if hipp==0
    % Load experiment table .........
    load ~/TT/HGexp.mat;
elseif hipp==1
    % if Hippo stimulation es-TT........
    load('~/TT/Hippoexp.mat');Hippoexp= [Hippoexp Hippoexp(:,3) Hippoexp(:,3)];
    HGexp=Hippoexp;
end

% select runs to average 
seltr = [7 8];
clear erbp; 
for m=1:length(seltr)
    n = seltr(m)
    pid = HGexp(n,1);  
    if HGexp(n,2)<100
            ss=['0',num2str(HGexp(n,2))];
    else
            ss= num2str(HGexp(n,2));
    end  
    side = HGexp(n,5);
    site = HGexp(n,4);

    disp(' ~~~~~~~~~~~~~~~~~ ' )
    fprintf('%s%1.0f%s%s\n','      Patient ',pid,' - Run ',ss);
    ss = ['~/TT/erbp/induce_estt_',num2str(pid),'_',ss,'.mat'];
    load(ss);
    if m==1
        data = erbp;
    else
        data = data+erbp;
    end
end
baseline = find(t>=-0.5 & t<=-0.2);
erbp = data/m;

% Apply threshold....
r = -0.55 + (-0.075+0.5)*rand(1000,1);
for k = 1:1000
    d = find(t>=r(k) & t<r(k)+0.005);
    basem(k,:) =mean(erbp(d,:),1);
end
TH= prctile(basem,[2.5 97.5],1);
erbpth = erbp;
for k=1:size(erbp,2)
    i = find(erbp(:,k)>=TH(1,k) & erbp(:,k)<=TH(2,k));
    erbpth(i,k)=0;
end

% load mesh ................
addpath /home/hiroyuki/gifti-1.8/;
addpath ~/afni_matlab/matlab
if side==1
    g= gifti('~/TT/icbm152/SUMA/std.60.lh.pial.gii');
    Bl=importdata('/home/hiroyuki/TT/icbm152/SUMA/std.60.lh.curv.niml.dset');
elseif side==0
    g= gifti('~/TT/icbm152/SUMA/std.60.rh.pial.gii');
    gp= gifti('~/TT/icbm152/SUMA/std.60.rh.sphere.gii');
    Bl=importdata('/home/hiroyuki/TT/icbm152/SUMA/std.60.rh.curv.niml.dset');
end
B = Bl.data;

%////////////////////////////////////////////////////////// %
% Draw Brain Surface 
addpath ~/es-fMRI
close all;figure;
renderFSmeshr_es(g.vertices, g.faces+0,B,[-90 0]);close(figure(1));
colormap(gray(1026));caxis([-0.8 2]);H=findobj(gcf,'type','patch');   
set(H,'edgealpha',0,'specularexponent',0.35,'specularstrength',0.2,...
'specularColorReflectance',1,'AmbientStrength',1,'Diffusestrength',0.5);
set(gcf,'color','w'); set(gca,'xcolor','w','ycolor','w','fontsize',24)
if side ==0 
    view(90,0);
    o=findobj(gca,'type','light');colormap(gray(1024));delete(o)
    light('position',[1000 0 0],'color','w');
    light('position',[100 150 0],'color','w');
    light('position',[100 -150 0],'color','w');
else
    view(-90,0);
    o=findobj(gca,'type','light');colormap(gray(1024));delete(o)
    light('position',[-1000 0 0],'color','w');
    light('position',[-100 150 0],'color','w');
     light('position',[-100 -150 0],'color','w');
end
% ---------------------------------------------------------  %
a=colorcet('CBD2');
colscale=a;
baseline = find(t>=-0.5 & t<=-0.2);

if hipp==0
    if site ==2
        ss = ['~/TT/erbp/AEP-Ave-ALHG_movie_estt_',num2str(pid),'.avi'];
    elseif site==1
        ss = ['~/TT/erbp/AEP-Ave-PMHG_movie_estt_',num2str(pid),'.avi'];
    elseif site ==3
        ss = ['~/TT/erbp/AEP-Ave-PlanumT_movie_estt_',num2str(pid),'.avi'];
    end    
elseif hipp==1
        ss = ['~/TT/erbp/AEP-Ave-Hippo_movie_estt_',num2str(pid),'.avi'];
end

v = VideoWriter(ss,'Motion JPEG AVI');
v.FrameRate = 8; 
v.Quality = 50;
maxz = 2.5;   % max/min colorscale value....
vv = floor(length(colscale)/(2*maxz)); 
tperiod = -0.0075:0.005:0.5;

open(v)
for ii = 1:length(tperiod)-1
    dur = find(t>=tperiod(ii) & t<tperiod(ii+1));
    ctime = floor(1000*(tperiod(ii)+tperiod(ii+1))/2);         

    if ctime < 0.004*1000
        x = zeros(length(CH),1);
    else
        x = mean(erbpth(dur,:),1);
    end
    title([num2str(ctime) ' ms'],'fontsize',36);
    o=findobj(gca,'type','line');delete(o);

    for k=1:size(GG,1)
        if x(k)>=0
            cc = floor(min(x(k)*vv,128))+128;
        else
            cc = floor(max(x(k)*vv,-128))+129;
        end
        if side==0
            h2=plot3(GG(k,1)+20, GG(k,2),GG(k,3),'o','markersize',28,...
                'markerfacecolor',[colscale(cc,:)],'markeredgecolor',[colscale(cc,:)]); 
            hold on; 
        else 
            h2=plot3(GG(k,1)-20, GG(k,2),GG(k,3),'o','markersize',28,...
                'markerfacecolor',[colscale(cc,:)],'markeredgecolor',[colscale(cc,:)]); 
            hold on; 
        end
    end      

    % Set colorbar property .....
    if ii==1
        yy=[];
        o=get(gcf,'ColorMap'); nl = size(o,1);
        yy(1:nl,:)=o; yy(nl+1:nl+256,:)=colscale;
        set(gcf,'colormap',yy); ca = caxis;H=linspace(ca(1), ca(2), size(yy,1));
        h=colorbar;  
        set(h,'Location','south','position',[0.4 0.15 0.22 0.02],'Limits',[H(1025) H(end)],...
            'ticks',[H(1025) H(end)],'ticklabels',[-maxz maxz],'fontsize',36);
    end
    drawnow;
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

%%



















