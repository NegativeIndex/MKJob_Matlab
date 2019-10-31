% clc
% clearvars
% close all

function gd_function(spara)

% wl=3.8;
% a=8;  % period
% b=1; % contact width
% t_etch=1.0;
% t_coat=0.5;
% 
% theta=0;
% phi=0;
% 
% sfolder='RawData';
% sig=sprintf('G2D_wl%0.2f_A%0.3fB%0.3f_TE%0.2fC%0.2f_Th%0.1fPh%0.1f',...
%     wl,a,b,t_etch,t_coat,theta,phi);
% spara=build_para_struct(wl,a,b,t_etch,t_coat,theta,phi,sig,sfolder);
% 
% 
% switch computer()
%     case 'PCWIN64'
%         rfolder='H:\functions';
%     case 'GLNXA64'
%         rfolder='/Users/wdai11/function_comsol';
%     otherwise
%         waring('Wrong comupter type')
% end
% addpath(fullfile(rfolder,'GD_Calc'));
% addpath(fullfile(rfolder,'Material_20190929'));
% addpath(fullfile(rfolder,'Matlab_20190607'));

% check the file existance
to_calculate=true;
to_display=(~to_calculate);  
  
fprintf(1,'%s to run\n',spara.sig);
if ~save_RTA_results(fullfile(spara.sfolder,spara.sig),0,0) && to_calculate
    return 
end



%% Define parameters for grating structure and incident field:
nm=1e-9;
um=1e-6;

wavelength=spara.wl;
theta=spara.theta*pi/180; % incidence polar angle from x1 axis
phi=spara.phi*pi/180; % incidence azimuthal angle around x1 axis (zero on x2 axis)
t_contact=0.1;
t_protect=0.1;
t_mirror=0.2;
mat_contact=spara.mat_contact;
mat_mirror=spara.mat_mirror;

% mat structure array
mat_set={'air','GaSb','Spacer','Al','Pd'};
nn_set=zeros(size(mat_set));
for i=1:length(mat_set)
    switch mat_set{i}
        case 'air'
            n=1;
        case 'GaSb'
            n=3.8;
        case 'Spacer'
            n=1.4;
        case {'Al','Pd'}
            n=getN_wl(sprintf('%s-Meep',mat_set{i}),wavelength*um);
        otherwise
            error('Material not defined')
    end
    nn_set(i)=n;        
end
mat2idx_map=containers.Map(mat_set,1:length(mat_set));
idx_sup=mat2idx_map('GaSb');
idx_spacer=mat2idx_map('Spacer');
idx_sub=mat2idx_map('air');


%%
m_max=50; % maximum diffraction order index

%% Construct grating.
clear grating
grating.pmt=cellfun(@(x)x^2, num2cell(nn_set),...
    'uniformoutput',false); % grating material permittivities
grating.pmt_sub_index=idx_sub; % substrate permittivity index
grating.pmt_sup_index=idx_sup; % superstrate permittivity index

% Define the x2 and x3 projections of the first grating period
% (d21,d31) and second grating period (d22,d32). The second period is
% parallel to the x2 axis.
grating.d21=spara.a; % first grating period: x2 projection
grating.d31=0;       % first grating period: x3 projection
grating.d22=0;       % second grating period: x2 projection
grating.d32=spara.a; % second grating period: x3 projection
grating.stratum={};

clear stratum stripe

% contour before coating without photolithography
x_contact=spara.b/2;
ptx1=[-spara.a/2,-x_contact,-x_contact];
pty1=[spara.t_etch,spara.t_etch,-(t_contact+t_protect)];
ptx1=[ptx1, -1*flip(ptx1)];
pty1=[pty1,flip(pty1)];

% spacer contour
if spara.t_coat>0.01
    [ptx2,pty2]=grow_downword(ptx1,pty1,spara.t_coat,1.5,16);
end
% mirror contour
[ptx3,pty3]=grow_downword(ptx1,pty1,spara.t_coat+t_mirror,1.5,16);
% plot contour
if to_display
    figure('position', [2000,450,560,420])
    hold on
    plot(ptx1,pty1,'-o')
    if spara.t_coat>0.01
        plot(ptx2,pty2,'-o')
    end
    plot(ptx3,pty3,'-o')
    hold off
    axis equal
    box on
end
% from bottom to top
if spara.t_coat>0.01
    bnds={{ptx3,pty3},{ptx2,pty2},{ptx1,pty1}};
    mats={'air',mat_mirror,'Spacer',mat_contact};
    mat_indices=values(mat2idx_map,mats);
    % polygons for air, Al mirror,spacer
    polygons=build_polygons_from_boundary(bnds,mat_indices,0.2);
else
    bnds={{ptx3,pty3},{ptx1,pty1}};
    mats={'air',mat_mirror,mat_contact};
    mat_indices=values(mat2idx_map,mats);
    % polygons for air, Al mirror
    polygons=build_polygons_from_boundary(bnds,mat_indices,0.2);
end
polygons(end)=[];

% contact polygon
polygons(end+1)=struct('polygon', ...
    polyshape(x_contact*[-1,1,1,-1],...
    [-t_contact,-t_contact,0,0]), ...
    'index',mat2idx_map(mat_contact));
% pretect polygon
polygons(end+1)=struct('polygon', ...
    polyshape(x_contact*[-1,1,1,-1],...
    -t_contact-[t_protect,t_protect,0,0]), ...
    'index',mat2idx_map(mat_mirror));
% GaSb pillar polygon
if spara.t_etch>0.01
    polygons(end+1)=struct('polygon', ...
        polyshape(x_contact*[-1,1,1,-1],...
        [spara.t_etch,spara.t_etch,0,0]), ...
        'index',mat2idx_map('GaSb'));
end
% GaSb substrate polygon
polygons(end+1)=struct('polygon', ...
    polyshape(spara.a/2*[-1,1,1,-1],...
    [spara.t_etch,spara.t_etch,spara.t_etch+0.2,spara.t_etch+0.2]), ...
    'index',mat2idx_map('GaSb'));
% plot polygons
if to_display
    figure('position', [2600,450,560,420])
    hold on
    for i=1:length(polygons)
        color=mat2color(mat_set{polygons(i).index});
        if isempty(color)
           color=[1,1,1]; 
        end
        % disp(color);
        plot(polygons(i).polygon,'facecolor',color);
    end
    hold off
    axis equal
    axis off
    box on
    figname=sprintf('Geom_A%0.3fB%0.3f_TE%0.2fC%0.2f.jpg',...
        spara.a,spara.b,spara.t_etch,spara.t_coat);
    print(figname,'-djpeg','-r600')
    disp(figname);
end
% build stratums
stratums=stratums_from_polygons(polygons);
grating.stratum=stratums;
clear stratum stripe

%% Define the indicent field.
clear inc_field
inc_field.wavelength=wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1=-cos(theta)/wavelength*nn.
inc_field.f2=sin(theta)*cos(phi)/wavelength*nn_set(idx_sup);
inc_field.f3=sin(theta)*sin(phi)/wavelength*nn_set(idx_sup);

% Specify which diffracted orders are to be retained in the calculations.
% The orders are specified so that the retained diffracted orders' spatial
% frequencies form a hexagonal pattern such as that illustrated in
% GD-Calc.pdf, equations 4.13 and 4.14 and Figure 5.
clear order
order(1).m2=0;
order(1).m1=-m_max:m_max;

%% Run the diffraction calculations.
if to_calculate
    
    [param_size,scat_field,inc_field]=gdc(grating,inc_field,order,false);
    
    % Compute the diffraction efficiencies.
    [R,T]=gdc_eff(scat_field,inc_field);
    % Discard diffracted waves that decay exponentially with distance from the
    % grating. (These include evanescent waves and, if the substrate's
    % permittivity is not real-valued, all transmitted waves.)
    R=R(imag([scat_field.f1r])==0);
    T=T(imag([scat_field.f1t])==0);
    % Tabulate the diffraction order indices and diffraction efficiencies for
    % an incident field polarized parallel to the incident plane. Also tabulate
    % the fractional energy loss in the grating.
    disp(' ');
    disp('Diffraction efficiencies (m1, m2, eff1, eff2, eff3, eff4)');
    disp('R:');
    disp(num2str([[R.m1].' [R.m2].' [R.eff1].' [R.eff2].' [R.eff3].' [R.eff4].']));
    disp('T:');
    disp(num2str([[T.m1].' [T.m2].' [T.eff1].' [T.eff2].' [T.eff3].' [T.eff3].']));
    disp('Energy loss:');
    disp([ num2str(1-sum([[[R.eff1].']; [[T.eff1].']])) ...
        , '   ', num2str(1-sum([[[R.eff2].']; [[T.eff2].']])) ...
        , '   ', num2str(1-sum([[[R.eff3].']; [[T.eff3].']])) ...
        , '   ', num2str(1-sum([[[R.eff4].']; [[T.eff4].']])) ...
        ]);
    
    save_RTA_results(fullfile(spara.sfolder,spara.sig),R,T);
end
%% Plot the grating. (Caution: Cancel the plot if N is very large.)

if to_display
    clear pmt_display
    gdc(grating);
    for i=1:length(mat_set)
        color=mat2color(mat_set{i});
        pmt_display(i).color=color;
        pmt_display(i).alpha=1;
        pmt_display(i).name=mat_set{i};
    end   
    x_limit=[-1.5*spara.a,-1.5*spara.a,-1.5*spara.a; ...
        1.5*spara.a,1.5*spara.a,1.5*spara.a];
    gdc_plot(grating,1,pmt_display,x_limit);
    view(0,4);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% end of main funciton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
function color=mat2color(mat)
switch mat
    case 'air'
        % color=[1,1,1];
        color=[];
    case 'Spacer'
        color=[0,1,0.5];
    case 'GaSb'
        color=[0.25,0.25,0.25];
    case {'Al'}
        color=[218,165,32]/255;
    case {'Pd'}
        color=[210,105,45]/255;
    otherwise
        error('Material not defined');
end

end

%%
function newcalc=save_RTA_results(sig,R,T)
% newcalc=true means we need run the simulations
% newcalc=false means all the files already exist and we don't do anything
fname_R=[sig,'_Refl.txt'];
fname_T=[sig,'_Tran.txt'];
fname_A=[sig,'_Absp.txt'];

if ~isstruct(R)  % if R is not struct, then check the existance of the three files
    newcalc=false;
    if ~exist(fname_R, 'file') newcalc=true; end
    if ~exist(fname_T, 'file') newcalc=true; end
    if ~exist(fname_A, 'file') newcalc=true; end
    if ~newcalc
        fprintf(1,'All files exist\n');
        return
    end
    
else  % save R,T,A
    fname=fname_R;
    fid=fopen(fname,'w');
    ss='%%%%Refl m1, m2, eff1, eff2 eff3 eff4\n';
    fprintf(fid,ss);
    ss=num2str([[R.m1].' [R.m2].' [R.eff1].' [R.eff2].' [R.eff3].' [R.eff4].']);
    for i=1:size(ss,1)
        fprintf(fid,'%s\n',ss(i,:));
    end
    fclose(fid);
    
    fname=fname_T;
    fid=fopen(fname,'w');
    ss='%%%%Tran: m1, m2, eff1, eff2 eff3 eff4\n';
    fprintf(fid,ss);
    ss=num2str([[T.m1].' [T.m2].' [T.eff1].' [T.eff2].' [T.eff3].' [T.eff4].']);
    for i=1:size(ss,1)
        fprintf(fid,'%s\n',ss(i,:));
    end
    fclose(fid);
    
    fname=fname_A;
    fid=fopen(fname,'w');
    ss='%%%%Absp: eff1, eff2 eff3 eff4\n';
    fprintf(fid,ss);
    ss1=num2str(1-sum([[[R.eff1].']; [[T.eff1].']]));
    ss2=num2str(1-sum([[[R.eff2].']; [[T.eff2].']]));
    ss3=num2str(1-sum([[[R.eff3].']; [[T.eff3].']]));
    ss4=num2str(1-sum([[[R.eff4].']; [[T.eff4].']]));
    fprintf(fid,'%s %s %s %s\n',ss1,ss2,ss3,ss4);
    fclose(fid);   
    fprintf(1,'Save all files\n');
    newcalc=true;
end

end


%% build polygon from material boundaris
function polygons_struct=build_polygons_from_boundary(bnds,mat_indices,ext)
% bnds is one item shorter than mat_indices
% bnds{i} is the boundary between material mat_indices(i) and
% mat_indices(i+1)
% polygon_struct include polygon and mat_indices
% polygons_struct is an array of strucutre
% ext is the region before the first boundary and after the last boundary

n=length(mat_indices);
polygons_struct=struct('polygon',cell(1,n),'index',cell(1,n));

% get boundary points
xx=[]; yy=[];
for i=1:length(bnds)
    xx=[xx,bnds{i}{1}];
    yy=[yy,bnds{i}{2}];
end

xmin=min(xx); xmax=max(xx);
ymin=min(yy); ymax=max(yy);

bnd2s=cell(1,length(bnds)+2);
bnd2s{1}={[xmin,xmax],[ymin-ext,ymin-ext]};
bnd2s(2:end-1)=bnds;
bnd2s{end}={[xmin,xmax],[ymax+ext,ymax+ext]};

for i=1:n
    % fprintf(1,'%d polygon\n',i);
    x1=bnd2s{i}{1}; x2=bnd2s{i+1}{1};
    y1=bnd2s{i}{2}; y2=bnd2s{i+1}{2};
    polygons_struct(i).polygon=polyshape([x1,flip(x2)],[y1,flip(y2)]);
    polygons_struct(i).index=mat_indices{i};
end

end

%% 
function stratums=stratums_from_polygons(polygons)

% get vertices 
xx=[]; yy=[];
for i=1:length(polygons)
    [xx0,yy0]=boundary(polygons(i).polygon);
    xx=[xx;xx0];yy=[yy; yy0];
end
xx=xx(~isnan(xx)); yy=yy(~isnan(yy));

xmin=min(xx); xmax=max(xx);
ymin=min(yy); ymax=max(yy);

tola=max(abs(yy));
tol=0.01;
yyticks=sort(uniquetol(yy,tol));
yyticks=fine_ticks(yyticks,0.3);
% cut stratums
stratums=cell(1,length(yyticks)-1);
clear stratum stripe
% stratum template
stratum.type=1; % uniperiodic stratum
% The following h11, h12 spec indicates that the stratum's period
% vector matches the first grating period (GD-Calc.pdf, equations 3.22
% and 3.23).
stratum.h11=1;
stratum.h12=0;

for i=1:length(stratums)
    y1=yyticks(i);
    y2=yyticks(i+1);
    stratum.thick=y2-y1; % stratum thickness
    % cut stripe
    % fprintf(1,'y1=%0.4f,y2=%0.4f\n ',y1,y2);
    
    rects_not_ready=true;
    repeat=0;
    while rects_not_ready 
        ymid=(y1+y2)/2+(rand()-0.5)*tol*(y2-y1)/10;
        g0=tol*(y2-y1)/100;
        midline=polyshape([xmin,xmax,xmax,xmin],...
            [ymid-g0,ymid-g0,ymid+g0,ymid+g0]);
        
        rects={};  % geometry
        rm_indices={};  % material
        for j=1:length(polygons)
            polyout=intersect(polygons(j).polygon,midline);
            rect0=regions(polyout);
            rmat0=cell(size(rect0));
            rmat0(:)={polygons(j).index};
            
            rects=[rects;rect0];
            rm_indices=[rm_indices;rmat0];
        end
        idx=sort_idx_rects(rects);
        rects=rects(idx);
        rm_indices=rm_indices(idx);
        % check rects should be monotonic
        rects_not_ready=is_overlap_rect(rects,1e-3);
        repeat=repeat+1;
        assert(repeat<10,'Geometry cannot be assembled properly');
    end
    
    % add stripe to stratum
    stratum.stripe={};
    for j=1:length(rects)
        stripe.pmt_index=rm_indices{j}; % first stripe's permittivity index
        [xleft,xright]=x_rect(rects(j));
        % fprintf(1,'[%0.4f,%0.4f] ',xleft,xright);
        stripe.c1=xright/(xmax-xmin); % first stripe's boundary on positive side
        stratum.stripe{j}=stripe;
    end
    % fprintf(1,'\n');
    stratums{i}=stratum;
end
clear stratum stripe

end

%% 
function out=fine_ticks(in,gap)
% add more points to array in so that the mesh is smaller than gap
    out=[];
    for i=1:length(in)-1
        out=[out,in(i)];
        n=ceil(abs(in(i+1)-in(i))/gap);
        if n>1
            temp=linspace(in(i),in(i+1),n+1);
            out=[out,temp(2:end-1)];
        end
        
    end
    out=[out,in(end)];
end

%% function about rects
function idx=sort_idx_rects(rects)
% sort rects based on centroid x
cenx=zeros(size(rects));
for k=1:length(rects)
    [cenx(k),~]=centroid(rects(k));
end
[~,idx]=sort(cenx);
end

function [xleft,xright]=x_rect(rect)
[x,~]=boundary(rect);
x=sort(x(1:end-1));
xleft=mean(x(1:2));
xright=mean(x(end-1:end));
end

function out=is_overlap_rect(rects,tol)

xx=zeros(1,2*length(rects));
for i=1:length(rects)
  [xx(2*i-1),xx(2*i)]=x_rect(rects(i));
end

out=false;
while length(xx)>2
   overlap=range_length(range_intersection(xx(1:2),xx(3:end)));
   if overlap>tol*range_length(xx)
       out=true;
       return
   end
   xx=xx(3:end);
end

end


%%
function [ptxx,ptyy]=grow_downword(ptx,pty,dd,kd,n)
% make sure ptx and pty can be cut into two symmetric parts
% then corner_grow_downword function
% dd is the deposition thickness
% kd is the deposition anisotropic, kd>1 means prefer vertical
% n is the number of lines
[xx1,yy1]=corner_grow_downword(ptx(1:3),pty(1:3),dd,kd,round(n/2));
[xx2,yy2]=corner_grow_downword(ptx(4:6),pty(4:6),dd,kd,round(n/2));
ptxx=[xx1,xx2];
ptyy=[yy1,yy2];

end


%%
function [ptxx,ptyy]=corner_grow_downword(ptx,pty,dd,kd,n)
% three point to define a corner, horizontal line is at the top
% dd is the deposition thickness
% kd is the deposition anisotropic
% n is the number of lines

d0=max(pty)-min(pty); % % origin thickness difference
xe=max(ptx)-min(ptx); % x extension
% deposition from the origin upwards
if dd>d0  % small arc
    yy=linspace(dd,dd-d0,n);
else      % quarter arc
    yy=linspace(dd,0,n);
end
xx=sqrt(dd^2-yy.^2)/kd;
if xx(end)>xe
    ye=sqrt(dd^2-(kd*xe).^2);
    yy=linspace(dd,ye,n);  % even smaller arc limited by x
    xx=sqrt(dd^2-yy.^2)/kd;
end

% assemble points
if abs(pty(3)-pty(2))> abs(pty(2)-pty(1))  % horizontal, then down
    ptxx=ptx(3)-xx(end:-1:1);
    ptyy=pty(3)-yy(end:-1:1);
    if abs(ptxx(1)-ptx(1))<xe/n % arc only
        ptxx(1)=ptx(1);
    elseif dd>d0*0.999  % no corner,add tolerence 
        ptxx=[ptx(1),ptxx];
        ptyy=[pty(1)-dd,ptyy];
    else  % with corner
        ptxx=[ptx(1),ptxx(1),ptxx];
        ptyy=[pty(1)-dd,pty(1)-dd,ptyy];
    end
else  % up, then horizontal
    % arc
    ptxx=ptx(1)+xx;
    ptyy=pty(1)-yy;
    if abs(ptxx(end)-ptx(3))<xe/n % arc only
        ptxx(end)=ptx(3);
    elseif dd>d0*0.999  % no corner
        ptxx=[ptxx,ptx(3)];
        ptyy=[ptyy,pty(3)-dd];
    else  % with corner
        ptxx=[ptxx,ptxx(end),ptx(3)];
        ptyy=[ptyy,pty(3)-dd,pty(3)-dd];
    end
end
% now the corner point is important
end
