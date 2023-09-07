clear all
close all

colourmap=repmat([90:2.2:255],[3,1])'./255;
nchan=68;
numpar=49;
%% Wrangle data
files=dir('/media/jdscasta/Jason_drive/Fingerprinting_2021/CAMCAN_oct2021/sub_*fooof_training.mat');
files02=dir('/media/jdscasta/Jason_drive/Fingerprinting_2021/CAMCAN_oct2021/sub_*fooof_validation.mat');

count=1;
for i=1:numpar
    
    load(strcat(files(i,:).folder,'/',files(i,:).name));
    Session1_Gaussian{count,1}=Options.FOOOF;
    
    load(strcat(files02(i,:).folder,'/',files02(i,:).name));
    Session2_Gaussian{count,1}=Options.FOOOF;

    count=count+1;
    
end 

channel=[RowNames'];


%% Plot PSD

files=dir('/media/jdscasta/Jason_drive/Fingerprinting_2021/CAMCAN_oct2021/sub_*_training.mat');
files02=dir('/media/jdscasta/Jason_drive/Fingerprinting_2021/CAMCAN_oct2021/sub_*_validation.mat');

files=files(2:2:numpar*2);
files02=files02(2:2:numpar*2);

count=1;
for i=1:numpar
    
    load(strcat(files(i,:).folder,'/',files(i,:).name));
    PSD_Session1(count,:,:)=squeeze(TF(:,:,1:301));
    
    load(strcat(files02(i,:).folder,'/',files02(i,:).name));
    PSD_Session2(count,:,:)=squeeze(TF(:,:,1:301));
    count=count+1;

end 


PSD_data.PSD_Session1=PSD_Session1;
PSD_data.PSD_Session2=PSD_Session2;
PSD_data.ROI=channel;
PSD_data.Freq=Freqs;

figure
clear g

g(1,1)=gramm('x',Freqs(2:301),'y', squeeze(mean(log(PSD_Session1(:,:,2:301)),1)), 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('PSD Session 1');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();

figure
clear g

g(1,1)=gramm('x',log(Freqs(2:301)),'y', squeeze(mean(log(PSD_Session1(:,:,2:301)),1)), 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('PSD Session 1');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();

figure
clear g

g(1,1)=gramm('x',Freqs(2:151),'y', squeeze(mean(log(PSD_Session2(:,:,2:151)),1)), 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('PSD Session 2');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();

figure
clear g

g(1,1)=gramm('x',log(Freqs(2:151)),'y', squeeze(mean(log(PSD_Session2(:,:,2:151)),1)), 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('PSD Session 2');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();



base_psd=squeeze(mean(mean(log(PSD_Session1(:,:,2:151)),1),2));
ind_psd=squeeze(mean(mean(log(PSD_Session2(:,:,2:151)),1),2));
figure
clear g

g(1,1)=gramm('x',repmat(Freqs(2:151),[2,1]),'y', [base_psd,ind_psd]', 'color',{'Session 1', 'Session 2'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('PSD across all electrodes');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();

%% MODEL FITS 
for i=1:numpar
   
    mse_Session1_gaussian(i,:)=  [Session1_Gaussian{i,1}.stats.MSE];
    rsq_Session1_gaussian(i,:)=  [Session1_Gaussian{i,1}.stats.r_squared];

     mse_Session2_gaussian(i,:)=  [Session2_Gaussian{i,1}.stats.MSE];
    rsq_Session2_gaussian(i,:)=  [Session2_Gaussian{i,1}.stats.r_squared];

end


% mean fit across models and weight vs no weight 
mse_Session1=[mean(mean(mse_Session1_gaussian))]
rsq_Session1=[mean(mean(rsq_Session1_gaussian))]

% plot and descriptive stats over fits 
mse_Session1=[mean(mse_Session1_gaussian,2)];
rsq_Session1=[mean(rsq_Session1_gaussian,2)];

mse_Session1=reshape(mse_Session1,[1,numpar]);
rsq_Session1=reshape(rsq_Session1,[1,numpar]);

c=cell(numpar,1);
c(:)={'Gaussian'};

labels={c{:}};

figure
clear g
g(1,1)=gramm('x',labels','y',mse_Session1,'color',labels');
%Raw data as scatter plot
%Boxplots
g(1,1).stat_boxplot();
g(1,1).set_title('stat_boxplot()');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','','y','MSE');
g.set_title('MSE across FOOOF methods');
g.axe_property('YLim',[0 0.1]);
%g.axe_property('XLim',[0 4]);
g.draw();

figure
clear g
g(1,1)=gramm('x',labels','y',rsq_Session1,'color',labels');
%Raw data as scatter plot
%Boxplots
g(1,1).stat_boxplot();
g(1,1).set_title('stat_boxplot()');
%These functions can be called on arrays of gramm objects
g(1,1).set_color_options('map','brewer2');
g.set_names('x','','y','Rsq');
g.set_title('R-Squared across FOOOF methods');
g.axe_property('YLim',[0.8 1]);
%g.axe_property('XLim',[0 4]);
g.draw();

data=[mean(rsq_Session1_gaussian,2);...
    mean(rsq_Session2_gaussian,2)];

c=cell(numpar,1);
c(:)={'Gaussian'};


labels={c{:},c{:}};

c=cell(numpar,1);
c(:)={'Session 1'};

c2=cell(numpar,1);
c2(:)={'Session 2'};
cond={c{:},c2{:}};

figure
clear g
g(1,1)=gramm('x',labels','y',data,'color',cond');
%Raw data as scatter plot
%Boxplots
g(1,1).stat_boxplot();
g(1,1).set_title('stat_boxplot()');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','','y','R Sq');
g.set_title('R sq');
g.axe_property('YLim',[0.7 1]);
%g.axe_property('XLim',[0 4]);
g.draw();



%% Frq wise Error 
for i=1:numpar
    
    temp1=Session1_Gaussian{i,1}.stats;
    temp4=Session2_Gaussian{i,1}.stats;

    
    for j=1:68
   
    fmse_Session1_gaussian(i,j,:)=  temp1(j).frequency_wise_error(1:78);
    
    fmse_Session2_gaussian(i,j,:)=  temp4(j).frequency_wise_error(1:78);
    
    
    end
end

figure
clear g

g(1,1)=gramm('x',repmat(Freqs(2:79),[68,1]),'y', squeeze(mean(fmse_Session1_gaussian,1)), 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Frequency-wise error Gaussian' );
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Error');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();


%% Foofed spectrum
for i=1:numpar
    
    temp1=Session1_Gaussian{i,1}.data;
    temp4=Session2_Gaussian{i,1}.data;


    for j=1:68
                
    fpsd_Session1_gaussian(i,j,:)=  temp1(j).fooofed_spectrum(1:78);
    
    fpsd_Session2_gaussian(i,j,:)=  temp4(j).fooofed_spectrum(1:78);
    
    
    end
end


gauss_psd=squeeze(mean(mean(log(fpsd_Session1_gaussian),1),2));

figure
clear g

g(1,1)=gramm('x',repmat(Freqs(2:79),[2,1]),'y', [base_psd(1:78,:), gauss_psd]', 'color',{'PSD','Gaussian'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('FOOOFed spectrum --Baseline');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();


gauss_psd=squeeze(mean(mean(log(fpsd_Session2_gaussian),1),2));
figure
clear g

g(1,1)=gramm('x',repmat(Freqs(2:79),[2,1]),'y', [base_psd(1:78,:), gauss_psd]', 'color',{'PSD','Gaussian'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('FOOOFed spectrum --Baseline');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();


gauss_psd=squeeze(mean(mean(log(fpsd_Session1_gaussian(:,:,:)),1),2));
Wgauss_psd=squeeze(mean(mean(log(fpsd_Session2_gaussian(:,:,:)),1),2));

figure
clear g

g(1,1)=gramm('x',repmat(Freqs(2:79),[2,1]),'y', [gauss_psd,Wgauss_psd]', 'color',{ 'S1 Gaussian', 'S2 Gaussian'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('PSD across all electrodes');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();


%% APERIODIC 
for i=1:numpar
   
    aper_Session1_gaussian{i,:}=  [Session1_Gaussian{i,1}.aperiodics];


    aper_Session2_gaussian{i,:}=  [Session2_Gaussian{i,1}.aperiodics];

   

end

%descriptive stats of aperiodic 
% See range of aperiodic values per channel 
count=1;
for i=1:numpar
    offest_Session1_gaussian(count,:)=[aper_Session1_gaussian{i,1}.offset];
    
    expon_Session1_gaussian(count,:)=[aper_Session1_gaussian{i,1}.exponent];
    
    offest_Session2_gaussian(count,:)=[aper_Session2_gaussian{i,1}.offset];
    
    expon_Session2_gaussian(count,:)=[aper_Session2_gaussian{i,1}.exponent];
    
    count=count+1;
    
end


% plot and descriptive stats over fits 

data=[mean(offest_Session1_gaussian,2);...
    mean(offest_Session2_gaussian,2)];

c=cell(numpar,1);
c(:)={'Gaussian'};


labels={c{:},c{:}};

c=cell(numpar,1);
c(:)={'Session 1'};

c2=cell(numpar,1);
c2(:)={'Session 2'};
cond={c{:},c2{:}};

figure
clear g
g(1,1)=gramm('x',labels','y',data,'color',cond');
%Raw data as scatter plot
%Boxplots
g(1,1).stat_boxplot();
g(1,1).set_title('stat_boxplot()');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','','y','Offset');
g.set_title('FOOOF Offset');
%g.axe_property('YLim',[0 0.1]);
%g.axe_property('XLim',[0 4]);
g.draw();


data=[mean(expon_Session1_gaussian,2);...
    mean(expon_Session2_gaussian,2)];

figure
clear g
g(1,1)=gramm('x',labels','y',data,'color',cond');
%Raw data as scatter plot
%Boxplots
g(1,1).stat_boxplot();
g(1,1).set_title('stat_boxplot()');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','','y','Exponent');
g.set_title('FOOOF Exponent');
%g.axe_property('YLim',[0 0.1]);
%g.axe_property('XLim',[0 4]);
g.draw();


%% alpha /beta

for i=1:numpar
   
    tempp=  [Session1_Gaussian{i,1}.peaks];
    tempp=struct2table(tempp);
    %alpha_Session1_gaussian{i,:,:}=tempp(find(strcmp(tempp{:,1}, 'Pz')),:)
    alpha_Session1_gaussian{i,:,:}=tempp(find(tempp{:,2} >= 8 &  tempp{:,2} <= 12),:);
    
    tempp=  [Session2_Gaussian{i,1}.peaks];
    tempp=struct2table(tempp);
    alpha_Session2_gaussian{i,:,:}=tempp(find(tempp{:,2} >= 8 &  tempp{:,2} <= 12),:);
    
end

for i=1:numpar
   
    tempp=  [Session1_Gaussian{i,1}.peaks];
    tempp=struct2table(tempp);
    %alpha_Session1_gaussian{i,:,:}=tempp(find(strcmp(tempp{:,1}, 'Pz')),:)
    beta_Session1_gaussian{i,:,:}=tempp(find(tempp{:,2} >= 15 &  tempp{:,2} <= 30),:);
    
    tempp=  [Session2_Gaussian{i,1}.peaks];
    tempp=struct2table(tempp);
    beta_Session2_gaussian{i,:,:}=tempp(find(tempp{:,2} >= 15 &  tempp{:,2} <= 30),:);
    
end

%% wrangle peak alpha into channel matrix 

count=1;
for i=1:numpar
    for j=1:nchan
        
    chanind=find(strcmp([alpha_Session1_gaussian{i,1}.channel], RowNames(j)));
    
    if ~ isnan(chanind)
        
        amplitude_alpha_Session1_gaussian(count,j)=[alpha_Session1_gaussian{i,1}(chanind,:).amplitude(1)];
        centralfreq_alpha_Session1_gaussian(count,j)=[alpha_Session1_gaussian{i,1}(chanind,:).center_frequency(1)];
        std_alpha_Session1_gaussian(count,j)=[alpha_Session1_gaussian{i,1}(chanind,:).std_dev(1)];

    else
       amplitude_alpha_Session1_gaussian(count,j)=NaN;
       centralfreq_alpha_Session1_gaussian(count,j)=NaN;
       std_alpha_Session1_gaussian(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
    chanind=find(strcmp([alpha_Session2_gaussian{i,1}.channel], RowNames(j)));
    
    if ~ isnan(chanind)
        
        amplitude_alpha_Session2_gaussian(count,j)=[alpha_Session2_gaussian{i,1}(chanind,:).amplitude(1)];
        centralfreq_alpha_Session2_gaussian(count,j)=[alpha_Session2_gaussian{i,1}(chanind,:).center_frequency(1)];
        std_alpha_Session2_gaussian(count,j)=[alpha_Session2_gaussian{i,1}(chanind,:).std_dev(1)];


    else
       amplitude_alpha_Session2_gaussian(count,j)=NaN;
       centralfreq_alpha_Session2_gaussian(count,j)=NaN;
       std_alpha_Session2_gaussian(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
    chanind=find(strcmp([beta_Session1_gaussian{i,1}.channel], RowNames(j)));
    
    if ~ isnan(chanind)
        
        amplitude_beta_Session1_gaussian(count,j)=[beta_Session1_gaussian{i,1}(chanind,:).amplitude(1)];
        centralfreq_beta_Session1_gaussian(count,j)=[beta_Session1_gaussian{i,1}(chanind,:).center_frequency(1)];
        std_beta_Session1_gaussian(count,j)=[beta_Session1_gaussian{i,1}(chanind,:).std_dev(1)];

    else
       amplitude_beta_Session1_gaussian(count,j)=NaN;
       centralfreq_beta_Session1_gaussian(count,j)=NaN;
       std_beta_Session1_gaussian(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
    chanind=find(strcmp([beta_Session2_gaussian{i,1}.channel], RowNames(j)));
    
    if ~ isnan(chanind)
        
        amplitude_beta_Session2_gaussian(count,j)=[beta_Session2_gaussian{i,1}(chanind,:).amplitude(1)];
        centralfreq_beta_Session2_gaussian(count,j)=[beta_Session2_gaussian{i,1}(chanind,:).center_frequency(1)];
        std_beta_Session2_gaussian(count,j)=[beta_Session2_gaussian{i,1}(chanind,:).std_dev(1)];


    else
       amplitude_beta_Session2_gaussian(count,j)=NaN;
       centralfreq_beta_Session2_gaussian(count,j)=NaN;
       std_beta_Session2_gaussian(count,j)=NaN;
        
    end
    end
    count=count+1;
end


%% Collect fits of aperiodic and periodic 


for i=1:numpar
    
    temp11=reshape([Session1_Gaussian{i,1}.data.ap_fit], [size([Session1_Gaussian{i,1}.data.ap_fit],2)/68,68]);

    fitaper_Session1_gaussian(i,:,:)=  temp11(1:78,:);
    
    temp11=reshape([Session1_Gaussian{i,1}.data.peak_fit], [size([Session1_Gaussian{i,1}.data.peak_fit],2)/68,68]);

    fitper_Session1_gaussian(i,:,:)=  temp11(1:78,:);

    temp11=reshape([Session2_Gaussian{i,1}.data.ap_fit], [size([Session2_Gaussian{i,1}.data.ap_fit],2)/68,68]);

    fitaper_Session2_gaussian(i,:,:)=  temp11(1:78,:);

    temp11=reshape([Session2_Gaussian{i,1}.data.peak_fit], [size([Session2_Gaussian{i,1}.data.peak_fit],2)/68,68]);

    fitper_Session2_gaussian(i,:,:)=  temp11(1:78,:);

end



%% Output



for i=1:numpar
    temp=squeeze(PSD_Session1(i,:,:));
    temp2=squeeze(PSD_Session2(i,:,:));
   
    csvwrite( strcat('CAMCAN_octspectra/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('CAMCAN_octspectra/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
    temp=squeeze(fitaper_Session1_gaussian(i,:,:))';
    temp2=squeeze(fitaper_Session2_gaussian(i,:,:))';
    
    csvwrite( strcat('CAMCAN_octaperiodic/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('CAMCAN_octaperiodic/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
    temp=squeeze(fitper_Session1_gaussian(i,:,:))';
    temp2=squeeze(fitper_Session2_gaussian(i,:,:))';
    
    csvwrite( strcat('CAMCAN_octperiodic/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('CAMCAN_octperiodic/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
    
    temp=squeeze(fmse_Session1_gaussian(i,:,:))';
    temp2=squeeze(fmse_Session2_gaussian(i,:,:))';
    
    csvwrite( strcat('CAMCAN_octresiduals/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('CAMCAN_octresiduals/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
end


%
% clear 
% close all
% 
% colourmap=repmat([90:2.2:255],[3,1])'./255;
% nchan=203;
% numpar=49;
% %% Wrangle data
% files=dir('/media/jdscasta/Jason_drive/CAMCAN_fingerprinting_newsenssor/sub_*fooof_training.mat');
% files02=dir('/media/jdscasta/Jason_drive/CAMCAN_fingerprinting_newsenssor/sub_*fooof_validation.mat');
% 
% load(strcat(files(1,:).folder,'/',files(1,:).name));
% grad_chan=readtable('CAMCAN_Channel.csv');
% grad_chan=table2cell(grad_chan);
% 
% [C,ia,ib] = intersect(RowNames,grad_chan, 'stable'); % find common chan
% 
% 
% count=1;
% for i=1:numpar
%     
%     load(strcat(files(i,:).folder,'/',files(i,:).name));
%     Session1_Gaussian{count,1}=Options.FOOOF;
%     
%     load(strcat(files02(i,:).folder,'/',files02(i,:).name));
%     Session2_Gaussian{count,1}=Options.FOOOF;
% 
%     count=count+1;
%     
% end 
% 
% channel=[grad_chan'];
% 
% 
% %% Plot PSD
% 
% files=dir('/media/jdscasta/Jason_drive/CAMCAN_fingerprinting_newsenssor/sub_*_training.mat');
% files02=dir('/media/jdscasta/Jason_drive/CAMCAN_fingerprinting_newsenssor/sub_*_validation.mat');
% 
% files=files(1:2:numpar*2);
% files02=files02(1:2:numpar*2);
% 
% 
% count=1;
% for i=1:numpar
%     
%     
%     load(strcat(files(i,:).folder,'/',files(i,:).name));
%     chan1(:,count)= [RowNames(ia)];
%     PSD_Session1(count,:,:)=squeeze(TF(ia,:,1:301));
%     
%     load(strcat(files02(i,:).folder,'/',files02(i,:).name));
%     PSD_Session2(count,:,:)=squeeze(TF(ia,:,1:301));
%     chan2(:,count)= [RowNames(ia)];
%     count=count+1;
% 
% end 
% 
% PSD_data.PSD_Session1=PSD_Session1;
% PSD_data.PSD_Session2=PSD_Session2;
% PSD_data.ROI=channel;
% PSD_data.Freq=Freqs;
% 
% figure
% clear g
% 
% g(1,1)=gramm('x',Freqs(2:151),'y', squeeze(mean(log(PSD_Session1(:,:,2:151)),1)), 'color', categorical(1:203));
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('PSD Session 1');
% g(1,1).no_legend();
% %g(1,1).set_color_options('map',colourmap);
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Log Power');
% %g.axe_property('YLim',[1e-15 1.5e-9]);
% %g.axe_property('XLim',[0 50]);
% g.draw();
% 
% figure
% clear g
% 
% g(1,1)=gramm('x',log(Freqs(2:151)),'y', squeeze(mean(log(PSD_Session1(:,:,2:151)),1)), 'color', categorical(1:203));
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('PSD Session 1');
% g(1,1).no_legend();
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Log Power');
% 
% g.draw();
% 
% figure
% clear g
% 
% g(1,1)=gramm('x',Freqs(2:151),'y', squeeze(mean(log(PSD_Session2(:,:,2:151)),1)), 'color', categorical(1:203));
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('PSD Session 2');
% g(1,1).no_legend();
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Log Power');
% g.draw();
% 
% figure
% clear g
% 
% g(1,1)=gramm('x',log(Freqs(2:151)),'y', squeeze(mean(log(PSD_Session2(:,:,2:151)),1)), 'color', categorical(1:203));
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('PSD Session 2');
% g(1,1).no_legend();
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Log Power');
% g.draw();
% 
% 
% 
% base_psd=squeeze(mean(mean(log(PSD_Session1(:,:,2:151)),1),2));
% ind_psd=squeeze(mean(mean(log(PSD_Session2(:,:,2:151)),1),2));
% figure
% clear g
% 
% g(1,1)=gramm('x',repmat(Freqs(2:151),[2,1]),'y', [base_psd,ind_psd]', 'color',{'Session 1', 'Session 2'});
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('PSD across all electrodes');
% g(1,1).set_color_options('map','brewer2');
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Log Power');
% %g.axe_property('YLim',[1e-15 1.5e-9]);
% %g.axe_property('XLim',[0 50]);
% g.draw();
% 
% %% MODEL FITS 
% for i=1:numpar
%    
%     mse_Session1_gaussian(i,:)=  [Session1_Gaussian{i,1}.stats.MSE];
%     rsq_Session1_gaussian(i,:)=  [Session1_Gaussian{i,1}.stats.r_squared];
% 
%     mse_Session2_gaussian(i,:)=  [Session2_Gaussian{i,1}.stats.MSE];
%     rsq_Session2_gaussian(i,:)=  [Session2_Gaussian{i,1}.stats.r_squared];
% 
% end
% 
% mse_Session1_gaussian = mse_Session1_gaussian(:,ia);
% rsq_Session1_gaussian = rsq_Session1_gaussian(:,ia);
% mse_Session2_gaussian = mse_Session2_gaussian(:,ia);
% rsq_Session2_gaussian = rsq_Session2_gaussian(:,ia);
% 
% % mean fit across models and weight vs no weight 
% mse_Session1=[mean(mean(mse_Session1_gaussian))]
% rsq_Session1=[mean(mean(rsq_Session1_gaussian))]
% 
% % plot and descriptive stats over fits 
% mse_Session1=[mean(mse_Session1_gaussian,2)];
% rsq_Session1=[mean(rsq_Session1_gaussian,2)];
% 
% mse_Session1=reshape(mse_Session1,[1,numpar]);
% rsq_Session1=reshape(rsq_Session1,[1,numpar]);
% 
% c=cell(numpar,1);
% c(:)={'Gaussian'};
% 
% labels={c{:}};
% 
% figure
% clear g
% g(1,1)=gramm('x',labels','y',mse_Session1,'color',labels');
% %Raw data as scatter plot
% %Boxplots
% g(1,1).stat_boxplot();
% g(1,1).set_title('stat_boxplot()');
% g(1,1).set_color_options('map','brewer2');
% %These functions can be called on arrays of gramm objects
% g.set_names('x','','y','MSE');
% g.set_title('MSE across FOOOF methods');
% g.axe_property('YLim',[0 0.1]);
% %g.axe_property('XLim',[0 4]);
% g.draw();
% 
% figure
% clear g
% g(1,1)=gramm('x',labels','y',rsq_Session1,'color',labels');
% %Raw data as scatter plot
% %Boxplots
% g(1,1).stat_boxplot();
% g(1,1).set_title('stat_boxplot()');
% %These functions can be called on arrays of gramm objects
% g(1,1).set_color_options('map','brewer2');
% g.set_names('x','','y','Rsq');
% g.set_title('R-Squared across FOOOF methods');
% g.axe_property('YLim',[0.8 1]);
% %g.axe_property('XLim',[0 4]);
% g.draw();
% 
% data=[mean(rsq_Session1_gaussian,2);...
%     mean(rsq_Session2_gaussian,2)];
% 
% c=cell(numpar,1);
% c(:)={'Gaussian'};
% 
% labels={c{:},c{:}};
% 
% c=cell(numpar,1);
% c(:)={'Session 1'};
% 
% c2=cell(numpar,1);
% c2(:)={'Session 2'};
% cond={c{:},c2{:}};
% 
% figure
% clear g
% g(1,1)=gramm('x',labels','y',data,'color',cond');
% %Raw data as scatter plot
% %Boxplots
% g(1,1).stat_boxplot();
% g(1,1).set_title('stat_boxplot()');
% g(1,1).set_color_options('map','brewer2');
% %These functions can be called on arrays of gramm objects
% g.set_names('x','','y','R Sq');
% g.set_title('R sq');
% g.axe_property('YLim',[0.7 1]);
% %g.axe_property('XLim',[0 4]);
% g.draw();
% 
% 
% 
% %% Frq wise Error 
% for i=1:numpar
%     
%     temp1=Session1_Gaussian{i,1}.stats;
%     temp4=Session2_Gaussian{i,1}.stats;
% 
%     
%     for j=1:306
%    
%     fmse_Session1_gaussian(i,j,:)=  temp1(j).frequency_wise_error(1:78);
%     
%     fmse_Session2_gaussian(i,j,:)=  temp4(j).frequency_wise_error(1:78);
%     
%     
%     end
% end
% 
% fmse_Session1_gaussian = fmse_Session1_gaussian(:,ia,:);
% fmse_Session2_gaussian = fmse_Session2_gaussian(:,ia,:);
% 
% figure
% clear g
% 
% g(1,1)=gramm('x',repmat(Freqs(2:79),[203,1]),'y', squeeze(mean(fmse_Session1_gaussian,1)), 'color', categorical(1:203));
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('Frequency-wise error Gaussian' );
% g(1,1).no_legend();
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Error');
% %g.axe_property('YLim',[1e-15 1.5e-9]);
% %g.axe_property('XLim',[0 50]);
% g.draw();
% 
% 
% %% Foofed spectrum
% for i=1:numpar
%     
%     temp1=Session1_Gaussian{i,1}.data;
%     temp4=Session2_Gaussian{i,1}.data;
% 
% 
%     for j=1:306
%                 
%     fpsd_Session1_gaussian(i,j,:)=  temp1(j).fooofed_spectrum(1:78);
%     
%     fpsd_Session2_gaussian(i,j,:)=  temp4(j).fooofed_spectrum(1:78);
%     
% 
%     end
% end
% 
% 
% fpsd_Session1_gaussian = fpsd_Session1_gaussian(:,ia,:);
% fpsd_Session2_gaussian = fpsd_Session2_gaussian(:,ia,:);
% 
% 
% gauss_psd=squeeze(mean(mean(log(fpsd_Session1_gaussian),1),2));
% 
% figure
% clear g
% 
% g(1,1)=gramm('x',repmat(Freqs(2:79),[2,1]),'y', [base_psd(1:78,:), gauss_psd]', 'color',{'PSD','Gaussian'});
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('FOOOFed spectrum --Baseline');
% g(1,1).set_color_options('map','brewer2');
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Log Power');
% %g.axe_property('YLim',[1e-15 1.5e-9]);
% %g.axe_property('XLim',[0 50]);
% g.draw();
% 
% 
% gauss_psd=squeeze(mean(mean(log(fpsd_Session2_gaussian),1),2));
% figure
% clear g
% 
% g(1,1)=gramm('x',repmat(Freqs(2:79),[2,1]),'y', [base_psd(1:78,:), gauss_psd]', 'color',{'PSD','Gaussian'});
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('FOOOFed spectrum --Baseline');
% g(1,1).set_color_options('map','brewer2');
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Log Power');
% %g.axe_property('YLim',[1e-15 1.5e-9]);
% %g.axe_property('XLim',[0 50]);
% g.draw();
% 
% 
% gauss_psd=squeeze(mean(mean(log(fpsd_Session1_gaussian(:,:,:)),1),2));
% Wgauss_psd=squeeze(mean(mean(log(fpsd_Session2_gaussian(:,:,:)),1),2));
% 
% figure
% clear g
% 
% g(1,1)=gramm('x',repmat(Freqs(2:79),[2,1]),'y', [gauss_psd,Wgauss_psd]', 'color',{ 'S1 Gaussian', 'S2 Gaussian'});
% %smooth plot 
% g(1,1).geom_line();
% g(1,1).set_title('PSD across all electrodes');
% g(1,1).set_color_options('map','brewer2');
% %These functions can be called on arrays of gramm objects
% g.set_names('x','Frequency (Hz)','y','Log Power');
% %g.axe_property('YLim',[1e-15 1.5e-9]);
% %g.axe_property('XLim',[0 50]);
% g.draw();
% 
% 
% %% APERIODIC 
% for i=1:numpar
%    
%     aper_Session1_gaussian{i,:}=  [Session1_Gaussian{i,1}.aperiodics];
% 
%     aper_Session2_gaussian{i,:}=  [Session2_Gaussian{i,1}.aperiodics];
% 
% end
% 
% %descriptive stats of aperiodic 
% % See range of aperiodic values per channel 
% count=1;
% for i=1:numpar
%     offest_Session1_gaussian(count,:)=[aper_Session1_gaussian{i,1}.offset];
%     
%     expon_Session1_gaussian(count,:)=[aper_Session1_gaussian{i,1}.exponent];
%     
%     offest_Session2_gaussian(count,:)=[aper_Session2_gaussian{i,1}.offset];
%     
%     expon_Session2_gaussian(count,:)=[aper_Session2_gaussian{i,1}.exponent];
%     
%     count=count+1;
%     
% end
% 
% offest_Session1_gaussian = offest_Session1_gaussian(:,ia);
% 
% expon_Session1_gaussian = expon_Session1_gaussian(:,ia);
% 
% offest_Session2_gaussian = offest_Session2_gaussian(:,ia);
% 
% expon_Session2_gaussian = expon_Session2_gaussian(:,ia);
% 
% % plot and descriptive stats over fits 
% 
% data=[mean(offest_Session1_gaussian,2);...
%     mean(offest_Session2_gaussian,2)];
% 
% c=cell(numpar,1);
% c(:)={'Gaussian'};
% 
% 
% labels={c{:},c{:}};
% 
% c=cell(numpar,1);
% c(:)={'Session 1'};
% 
% c2=cell(numpar,1);
% c2(:)={'Session 2'};
% cond={c{:},c2{:}};
% 
% figure
% clear g
% g(1,1)=gramm('x',labels','y',data,'color',cond');
% %Raw data as scatter plot
% %Boxplots
% g(1,1).stat_boxplot();
% g(1,1).set_title('stat_boxplot()');
% g(1,1).set_color_options('map','brewer2');
% %These functions can be called on arrays of gramm objects
% g.set_names('x','','y','Offset');
% g.set_title('FOOOF Offset');
% %g.axe_property('YLim',[0 0.1]);
% %g.axe_property('XLim',[0 4]);
% g.draw();
% 
% 
% data=[mean(expon_Session1_gaussian,2);...
%     mean(expon_Session2_gaussian,2)];
% 
% figure
% clear g
% g(1,1)=gramm('x',labels','y',data,'color',cond');
% %Raw data as scatter plot
% %Boxplots
% g(1,1).stat_boxplot();
% g(1,1).set_title('stat_boxplot()');
% g(1,1).set_color_options('map','brewer2');
% %These functions can be called on arrays of gramm objects
% g.set_names('x','','y','Exponent');
% g.set_title('FOOOF Exponent');
% %g.axe_property('YLim',[0 0.1]);
% %g.axe_property('XLim',[0 4]);
% g.draw();
% 
% 
% 
% %% Collect fits of aperiodic and periodic 
% 
% for i=1:numpar
%     
%     temp11=reshape([Session1_Gaussian{i,1}.data.ap_fit], [size([Session1_Gaussian{i,1}.data.ap_fit],2)/306,306]);
% 
%     fitaper_Session1_gaussian(i,:,:)=  temp11(1:78,:);
%     
%     temp11=reshape([Session1_Gaussian{i,1}.data.peak_fit], [size([Session1_Gaussian{i,1}.data.peak_fit],2)/306,306]);
% 
%     fitper_Session1_gaussian(i,:,:)=  temp11(1:78,:);
% 
%     temp11=reshape([Session2_Gaussian{i,1}.data.ap_fit], [size([Session2_Gaussian{i,1}.data.ap_fit],2)/306,306]);
% 
%     fitaper_Session2_gaussian(i,:,:)=  temp11(1:78,:);
% 
%     temp11=reshape([Session2_Gaussian{i,1}.data.peak_fit], [size([Session2_Gaussian{i,1}.data.peak_fit],2)/306,306]);
% 
%     fitper_Session2_gaussian(i,:,:)=  temp11(1:78,:);
% 
% end
% 
% fitper_Session1_gaussian = fitper_Session1_gaussian(:,:,ia);
% fitaper_Session1_gaussian = fitaper_Session1_gaussian(:,:,ia);
% 
% fitper_Session2_gaussian = fitper_Session2_gaussian(:,:,ia);
% fitaper_Session2_gaussian = fitaper_Session2_gaussian(:,:,ia);
% 
% 
% %% Output
% 
% 
% 
% for i=1:numpar
%     temp=squeeze(PSD_Session1(i,:,:));
%     temp2=squeeze(PSD_Session2(i,:,:));
%    
%     csvwrite( strcat('CAMCAN_newspectra_sensor/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
%     csvwrite(strcat('CAMCAN_newspectra_sensor/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
%     
%     temp=squeeze(fitaper_Session1_gaussian(i,:,:))';
%     temp2=squeeze(fitaper_Session2_gaussian(i,:,:))';
%     
%     csvwrite( strcat('CAMCAN_newaperiodic_sensor/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
%     csvwrite(strcat('CAMCAN_newaperiodic_sensor/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
%     
%     temp=squeeze(fitper_Session1_gaussian(i,:,:))';
%     temp2=squeeze(fitper_Session2_gaussian(i,:,:))';
%     
%     csvwrite( strcat('CAMCAN_newperiodic_sensor/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
%     csvwrite(strcat('CAMCAN_newperiodic_sensor/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
%     
%     
%     temp=squeeze(fmse_Session1_gaussian(i,:,:))';
%     temp2=squeeze(fmse_Session2_gaussian(i,:,:))';
%     
%     csvwrite( strcat('CAMCAN_newresiduals_sensor/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
%     csvwrite(strcat('CAMCAN_newresiduals_sensor/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
%     
% end
