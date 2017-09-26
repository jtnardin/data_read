%create_ind_cell_fret_data.m written 9-22-17 by JTN to take the individual cell
%and fret data cell from cell_profile_data.mat and compute and save the aligned
%individual cell profiles for the various situations (WT/shAcat) and
%(mock/EGF).

clear all; clc

load('cell_profile_data.mat')
load('fret_profile_data.mat')

%Now create individual cell profiles
ind_cell_data = cell(8,2);
ind_cell_data_sv = cell(8,2);

ind_fret_data = cell(8,2);
ind_fret_data_sv = cell(8,2);

for i = 1:8
    for j = 1:2
        
        %figure out minimum duration for triplicate data
        tend_mean = min(tend((j-1)*3+1:j*3,i));
        
        LE_init = zeros(3,1);
        
        %find three initial LE locations
        for k = 1:3
            LE_init(k) = leading_edge_calc(cell_data_1d_mod{(j-1)*3+k,i}(:,1),1:540,0.5,0);
        end
        
        %sort by initial LE location for alignment
        [LE_init,I]  =sort(LE_init);
        %align back of sheets to smallest initial LE
        cell_data_1d_mod{(j-1)*3+I(2),i}(1:LE_init(2) - LE_init(1) ,:) = [];
        cell_data_1d_mod{(j-1)*3+I(3),i}(1:LE_init(3) - LE_init(1) ,:) = [];
        %align front of sheets to largest initial LE
        cell_data_1d_mod{(j-1)*3+I(2),i}(end-(LE_init(3) - LE_init(2)) + 1:end,:) = [];
        cell_data_1d_mod{(j-1)*3+I(1),i}(end-(LE_init(3) - LE_init(1)) + 1:end,:) = [];
        
        %remove the same points from the FRET data
        fret_data_1d_mod{(j-1)*3+I(2),i}(1:LE_init(2) - LE_init(1) ,:) = [];
        fret_data_1d_mod{(j-1)*3+I(3),i}(1:LE_init(3) - LE_init(1) ,:) = [];
        %align front of sheets to largest initial LE
        fret_data_1d_mod{(j-1)*3+I(2),i}(end-(LE_init(3) - LE_init(2)) + 1:end,:) = [];
        fret_data_1d_mod{(j-1)*3+I(1),i}(end-(LE_init(3) - LE_init(1)) + 1:end,:) = [];
        
        %how many spatial points do we have now
        xend = size(cell_data_1d_mod{(j-1)*3+1,i},1);
                
        %initialize individual data entries
        ind_cell_data{i,j} = zeros(3,xend,tend_mean);
        ind_cell_data_sv{i,j} = zeros(xend,tend_mean);
        ind_fret_data{i,j} = zeros(3,xend,tend_mean);
        ind_fret_data_sv{i,j} = zeros(xend,tend_mean);
        
        %loop through triplicate data and each time point to smooth the
        %cell profile in space for each example. Then calculate
        %correspdonding std 
        for k = 1:tend_mean
            for l = 1:3
                ind_cell_data{i,j}(l,:,k) = smooth(cell_data_1d_mod{(j-1)*3+l,i}(:,k))';
                ind_fret_data{i,j}(l,:,k) = smooth(fret_data_1d_mod{(j-1)*3+l,i}(:,k))';
            end
            ind_cell_data_sv{i,j}(:,k) = std(ind_cell_data{i,j}(:,:,k),1);
            ind_fret_data_sv{i,j}(:,k) = std(ind_fret_data{i,j}(:,:,k),1);
        end
        
        
    end
end

%save
save('ind_cell_prof_data.mat','ind_cell_data','ind_cell_data_sv',...
    'ind_fret_data','ind_fret_data_sv')

