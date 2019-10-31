clc
clearvars
close all
% low loss contact
sfolder='RawData';
if ~exist(sfolder, 'dir')
    mkdir(sfolder)
end
% oldfiles=FileList(sfolder,'.*txt','withfolder');
% for i=1:length(oldfiles)
%     fprintf(1,'Delete %s\n',oldfiles{i});
%     delete(oldfiles{i});
% end

switch computer()
    case 'PCWIN64'
        rfolder='H:\functions';
    case 'GLNXA64'
        rfolder='/Users/wdai11/function_comsol';
    otherwise
        waring('Wrong comupter type')
end
addpath(fullfile(rfolder,'GD_Calc'));
addpath(fullfile(rfolder,'Material_20190929'));
addpath(fullfile(rfolder,'Matlab_20190607'));

%% mateirals
% GaSb: 3.8; AlGaSb: 3.5; Si: 3.4; SiO2: 1.4
% ZrO2: 2.05; MgF2: 1.35 

tstart=tic;
%% parameters
% unit is um
mat_contact='Pd';
mat_mirror='Al';
wl=3.8;
a=4;
b=0.5;
t_etch=0.6;
t_coat=1.0;

% end of parameters

% error('stop');
%% run a simulation for all angles
theta_list= [1:2:89];
phi_list= [0:2:90];
n=0;
for theta=theta_list
    for phi=phi_list
        n=n+1;
    end
end
fprintf(1,'There are %d simulations in this folder\n',n);
disp('========================');
for theta=theta_list
    for phi=phi_list
        sig=sprintf('G2D_C%s_M%s_wl%0.2f_A%0.3fB%0.3f_TE%0.2fC%0.2f_Th%0.1fPh%0.1f',...
            mat_contact,mat_mirror,wl,a,b,t_etch,t_coat,theta,phi);
        spara=build_para_struct(mat_contact,mat_mirror,...
            wl,a,b,t_etch,t_coat,theta,phi,sig,sfolder);
        
        fprintf(1,'%s\n',sig);
        t_one=tic;
        gd_function(spara);
        toc(t_one)
        toc(tstart)
        fprintf(1,'%s finished\n',sig);
        % error('stop')
        disp('========================');
    end
end


%% 
function permu=get_all_permutation(para_cells)
    index_cell=cell(size(para_cells));
    for i=1:length(index_cell)
       index_cell{i}=1:length(para_cells{i}); 
    end
    permu=combvec(index_cell{:});
end
