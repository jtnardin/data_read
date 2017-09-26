%read_in_cell_data.m written 9-22-17 by JTN to import and save data from
%experimental images and create 1d cell profiles. This data can then be 
%used to calculate mean profiles or individual profiles in separate files.


clear all; clc


%initialize data cell
cell_data_1d_mod = cell(6,8);

x = linspace(0,1,540);

%to determine regions of interest -- calculated by eye from experimental
%videos to see what regions are relevant, and for what time period.
tend = 144*ones(6,8);
xstart = 1*ones(6,8);
xend = 540*ones(6,8);
yend = 540*ones(6,8);

xstart(4,1) = 320;
xstart(4,2) = 311;
xstart(4,4) = 148;
xstart(5,2) = 230;
xstart(5,3) = 153;
xstart(6,2) = 37;
xstart(6,3) = 115;
xstart(6,4) = 39;

xend(4,2) = 505;
xend(4,4) = 503;

yend(4,2) = 500;
yend(4,3) = 472;
yend(5,4) = 490;

tend(4,1) = 125;
tend(4,3) = 112;
tend(5,3) = 131;
tend(5,4) = 104;
tend(6,2) = 123;
tend(6,3) = 112;
tend(6,4) = 112;

%loop through each video and get relevant data
for i = 1:6
    for j = 1:8
        %load in corresponding video
        fname = ['WL_2_' char(65+i) '0' num2str(j+1) '.tif'];

        if i == 5 && j == 3 %had to modify this video -- lots of noise in wound
            fname = ['WL_2_' char(65+i) '0' num2str(j+1) '_concat.tif'];
        end
            
        %retrive info of image
        info = imfinfo(fname);

        %number of images
        num_images = tend(i,j);

        %xnum, ynum correspond to height and width of video resp.
        xnum = info.Width;
        ynum = xend(i,j) - xstart(i,j) + 1;

        %initialize matrices for videos
        A = zeros(ynum,xnum,num_images);

        B = zeros(540,540);
                
        for k = 1:num_images
            B = imread(fname,k,'info',info);
            A(:,:,k) = B(xstart(i,j):xend(i,j),:);
        end
        
        %ignore errant cells in wound space
        if yend(i,j)~=540
            A(:,yend(i,j):540,:)=0; 
        end
                
        %normalizing factor gotten from data in back of sheet
        img_mean = squeeze(mean(mean(A(:,1:100,:))));

        
        %now put 1d info into big vector
        %initialize entry in data cell
        cell_data_1d_mod{i,j} = zeros(xnum,num_images);
        
        
        %remove any noise past the LE
        % do so by defining cutoff density
        if i == 5 && j ==3
            LE_dens = 0.06; %this sample is particularly noisy
        else
            LE_dens = 0.02;
        end

        for k = 1:num_images
            %normalized cell density
            cell_data_1d_mod{i,j}(:,k) = mean(A(:,:,k))/img_mean(k);
            %find where cell dens small
            LE = leading_edge_calc(cell_data_1d_mod{i,j}(:,k),x,LE_dens,0);
            %set those points to zero
            cell_data_1d_mod{i,j}(x>LE,k) = 0; 
        end
        
    end
end

% save('cell_profile_data.mat','cell_data_1d_mod','tend','xstart','yend','xend')
