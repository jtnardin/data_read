clear all; clc

load('ind_cell_prof_data.mat')

well = 3;

rep = 1;

for k = 1:144
    yyaxis left
    hold off
    plot(ind_cell_data{well,2}(rep,:,k),'b')
    yyaxis right
    plot(ind_fret_data{well,2}(rep,:,k),'r')
    hold on
    plot(ind_fret_data{well,2}(rep,:,k)./ind_cell_data{well,2}(rep,:,k),'color'...
        ,[0 .5 0])
    axis([0 size(ind_fret_data{well,2},2) 0 1.2*max(max(ind_fret_data{well,2}(1,:,:)))])
    pause(.125)
end