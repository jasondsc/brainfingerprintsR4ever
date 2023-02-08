%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%             Script to organize FOOOF & spectral
%       outputs from CAMCAN dataset
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

colourmap=repmat([90:2.2:255],[3,1])'./255;
nchan=68;
numpar=606;
%% Wrangle data
files_RS=dir('~/Documents/CAMCAN_outputs/CAMCAN_AUG_FOOOF/Rest*.mat');
files_Task=dir('~/Documents/CAMCAN_outputs/CAMCAN_AUG_FOOOF/Task*.mat');

count=1;
for i=1:numpar
    
    load(strcat(files_RS(i,:).folder,'/',files_RS(i,:).name));
    Rest{count,1}=Options.FOOOF;
    PSD_rest(count,:,:)=squeeze(TF);
    Ids{count,1}=extractBetween(files_RS(i,:).name, 1,12);
    
    load(strcat(files_Task(i,:).folder,'/',files_Task(i,:).name));
    Task{count,1}=Options.FOOOF;
    PSD_task(count,:,:)=squeeze(TF);
    %PSD_task(count,:,:)=squeeze(TF(:,:,1:8));

    count=count+1;
    
end 

aa=reshape(permute(PSD_rest(:,:,1:300), [1,3,2]), [606,68*300]);
sub1=squeeze(PSD_rest(1,:,1:300));

channel=[RowNames'];


PSD_data.PSD_rest=PSD_rest;
PSD_data.PSD_task=PSD_task;
PSD_data.ROI=channel;
PSD_data.Freq=Freqs;


figure
clear g

g(1,1)=gramm('x',Freqs(2:301),'y', squeeze(mean(log(PSD_rest(:,:,2:301)),1)), 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('PSD Rest');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_full_specra_rest.pdf')


figure
clear g

g(1,1)=gramm('x',log(Freqs(2:301)),'y', squeeze(mean(log(PSD_rest(:,:,2:301)),1)), 'color', categorical(1:68));
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

g(1,1)=gramm('x',Freqs(2:301),'y', squeeze(mean(log(PSD_task(:,:,2:301)),1)), 'color', categorical(1:68));
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
saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_full_specra_task.pdf')


figure
clear g

g(1,1)=gramm('x',log(Freqs(2:301)),'y', squeeze(mean(log(PSD_task(:,:,2:301)),1)), 'color', categorical(1:68));
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
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/PSD_task_loglog.pdf')


base_psd=squeeze(mean(mean(log(PSD_rest(:,:,2:301)),1),2));
ind_psd=squeeze(mean(mean(log(PSD_task(:,:,2:301)),1),2));
figure
clear g

g(1,1)=gramm('x',repmat(Freqs(2:301),[2,1]),'y', [base_psd,ind_psd]', 'color',{'Rest', 'Task'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('PSD across all electrodes');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_full_specra_both.pdf')

%% MODEL FITS 
for i=1:numpar
   
    mse_Rest(i,:)=  [Rest{i,1}.stats.MSE];
    rsq_Rest(i,:)=  [Rest{i,1}.stats.r_squared];
    mse_Task(i,:)=  [Task{i,1}.stats.MSE];
    rsq_Task(i,:)=  [Task{i,1}.stats.r_squared];

end

figure
histogram(rsq_Task)
hold on;
histogram(rsq_Rest)

% plot and descriptive stats over fits 
mse=[mean(mse_Rest,2), mean(mse_Task,2)];
rsq=[mean(rsq_Rest,2),mean(rsq_Task,2)];

mse_Session1=reshape(mse,[1,2*numpar]);
rsq_Session1=reshape(rsq,[1,2*numpar]);

c=cell(numpar,1);
c(:)={'Rest'};

c2=cell(numpar,1);
c2(:)={'Task'};

labels={c{:},c2{:}};

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
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/MSE_Session1_methods.pdf')


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
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/RSQ_Session1_methods.pdf')


%% Frq wise Error 
for i=1:numpar
    
    temp1=Rest{i,1}.stats;
    temp2=Task{i,1}.stats;

    
    for j=1:68
   
    fmse_Rest(i,j,:)=  temp1(j).frequency_wise_error(1:78);
    fmse_Task(i,j,:)=  temp2(j).frequency_wise_error(1:78);

    
    end
end

figure
clear g

g(1,1)=gramm('x',repmat(Freqs(3:80),[68,1]),'y', squeeze(mean(fmse_Rest,1)), 'color', categorical(1:68));
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
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/fre_wise_error_gaussian.pdf')

figure
clear g

g(1,1)=gramm('x',repmat(Freqs(3:80),[68,1]),'y', squeeze(mean(fmse_Task,1)), 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Frequency-wise error Cauchy');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Error');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/fre_wise_error_cauchy.pdf')

figure
clear g

g(1,1)=gramm('x',repmat(Freqs(3:80),[68,1]),'y', squeeze(mean(fmse_Rest,1)), 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Frequency-wise error Mixed');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Error');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/fre_wise_error_mixed.pdf')

data=[squeeze(mean(mean(fmse_Rest,1))),squeeze(mean(mean(fmse_Task,1)))]';

figure
clear g

g(1,1)=gramm('x',repmat(Freqs(3:80),[2,1]),'y', data, 'color', {'Rest','Task'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Frequency-wise error across electrodes');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Error');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/fre_wise_error_compare.pdf')

%% Foofed spectrum
for i=1:numpar
    
    temp1=Rest{i,1}.data;
    temp2=Task{i,1}.data;

    for j=1:68
                
    fpsd_Rest(i,j,:)=  temp1(j).fooofed_spectrum(1:78);
    fpsd_Task(i,j,:)=  temp2(j).fooofed_spectrum(1:78);
    
    
    end
end


Rest_psd=squeeze(mean(mean(log(fpsd_Rest),1),2));
Task_psd=squeeze(mean(mean(log(fpsd_Task),1),2));

figure
clear g


g(1,1)=gramm('x',repmat(Freqs(3:80),[4,1]),'y', [squeeze(mean(mean(log(PSD_rest(:,:,3:80)),1),2)), Rest_psd, squeeze(mean(mean(log(PSD_task(:,:,3:80)),1),2)),Task_psd]', 'color',{'PSD Rest','FOOOF Rest', 'PSD Task', 'FOOOF Task'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('FOOOFed spectrum --Baseline');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_fooofed_specra_both.pdf')


% plot example good psd
figure
plot(Freqs(1:300), squeeze(log(PSD_rest(97,:,1:300)))')

% sub_CC110174
figure
plot(Freqs(1:300), squeeze(log(PSD_rest(20,:,1:300)))')

% sub_CC122405
figure
plot(Freqs(1:300), squeeze(log(PSD_rest(186,:,1:300)))')

%% APERIODIC 
for i=1:numpar
   
    aper_Rest{i,:}=  [Rest{i,1}.aperiodics];
    aper_Task{i,:}=  [Task{i,1}.aperiodics];
   
end

%descriptive stats of aperiodic 
% See range of aperiodic values per channel 
count=1;
for i=1:numpar
    offest_Rest(count,:)=[aper_Rest{i,1}.offset];
    offest_Task(count,:)=[aper_Task{i,1}.offset];
    
    expon_Rest(count,:)=[aper_Rest{i,1}.exponent];
    expon_Task(count,:)=[aper_Task{i,1}.exponent];

    count=count+1;
    
end


% plot and descriptive stats over fits 

data=[mean(offest_Rest,2);mean(offest_Task,2)];

c=cell(numpar,1);
c(:)={'Rest'};

c2=cell(numpar,1);
c2(:)={'Task'};

labels={c{:},c2{:}};

figure
clear g
g(1,1)=gramm('x',labels','y',data,'color',labels');
%Raw data as scatter plot
%Boxplots
g(1,1).stat_boxplot();
g(1,1).set_title('stat_boxplot()');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','','y','Offset');
g.set_title('FOOOF Offset');
g.axe_property('YLim',[-1.5 1.5]);
%g.axe_property('XLim',[0 4]);
g.draw();
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/fooofed_offset_compare.pdf')


data=[mean(expon_Rest,2);mean(expon_Task,2)];

figure
clear g
g(1,1)=gramm('x',labels','y',data,'color',labels');
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
%saveas(gcf,'/Users/jason/Desktop/FOOOF/Results/fooofed_expon_compare.pdf')


%% alpha /beta
for i=1:numpar

    tempp=  [Rest{i,1}.peaks];
   if length(tempp) ==0
        alpha_Rest{i,:,:}=[];
   else
    tempp=struct2table(tempp);
    %alpha_Session1_gaussian{i,:,:}=tempp(find(strcmp(tempp{:,1}, 'Pz')),:)
    alpha_Rest{i,:,:}=tempp(find(tempp{:,2} >= 8 &  tempp{:,2} <= 15),:);
   end
   
       tempp=  [Task{i,1}.peaks];
    if length(tempp) ==0
        alpha_Task{i,:,:}=[];
    else
    tempp=struct2table(tempp);
    alpha_Task{i,:,:}=tempp(find(tempp{:,2} >=8 &  tempp{:,2} <= 15),:);
    end
    
end
length(find(cellfun(@height, alpha_Rest)==0))
length(find(cellfun(@height, alpha_Task)==0))

for i=1:numpar
   
    tempp=  [Rest{i,1}.peaks];
   if length(tempp) ==0
        beta_Rest{i,:,:}=[];
   else
    tempp=struct2table(tempp);
    %alpha_Session1_gaussian{i,:,:}=tempp(find(strcmp(tempp{:,1}, 'Pz')),:)
    beta_Rest{i,:,:}=tempp(find(tempp{:,2} >= 15 &  tempp{:,2} < 30),:);
   end
   
       tempp=  [Task{i,1}.peaks];
    if length(tempp) ==0
        beta_Task{i,:,:}=[];
    else
    tempp=struct2table(tempp);
    beta_Task{i,:,:}=tempp(find(tempp{:,2} >= 15 &  tempp{:,2} < 30),:);
    end
        
end
length(find(cellfun(@height, beta_Rest)==0))
length(find(cellfun(@height, beta_Task)==0))


for i=1:numpar
   
    tempp=  [Rest{i,1}.peaks];
   if length(tempp) ==0
        delta_Rest{i,:,:}=[];
   else
    tempp=struct2table(tempp);
    %alpha_Session1_gaussian{i,:,:}=tempp(find(strcmp(tempp{:,1}, 'Pz')),:)
    delta_Rest{i,:,:}=tempp(find(tempp{:,2} >= 0 &  tempp{:,2} < 4),:);
   end
   
       tempp=  [Task{i,1}.peaks];
    if length(tempp) ==0
        delta_Task{i,:,:}=[];
    else
    tempp=struct2table(tempp);
    delta_Task{i,:,:}=tempp(find(tempp{:,2} >= 0 &  tempp{:,2} < 4),:);
    end
        
end

for i=1:numpar
   
    tempp=  [Rest{i,1}.peaks];
   if length(tempp) ==0
        theta_Rest{i,:,:}=[];
   else
    tempp=struct2table(tempp);
    %alpha_Session1_gaussian{i,:,:}=tempp(find(strcmp(tempp{:,1}, 'Pz')),:)
    theta_Rest{i,:,:}=tempp(find(tempp{:,2} >= 4 &  tempp{:,2} < 7),:);
   end
   
       tempp=  [Task{i,1}.peaks];
    if length(tempp) ==0
        theta_Task{i,:,:}=[];
    else
    tempp=struct2table(tempp);
    theta_Task{i,:,:}=tempp(find(tempp{:,2} >= 4 &  tempp{:,2} < 7),:);
    end
        
end

for i=1:numpar
   
    tempp=  [Rest{i,1}.peaks];
   if length(tempp) ==0
        gamma_Rest{i,:,:}=[];
   else
    tempp=struct2table(tempp);
    %alpha_Session1_gaussian{i,:,:}=tempp(find(strcmp(tempp{:,1}, 'Pz')),:)
    gamma_Rest{i,:,:}=tempp(find(tempp{:,2} >= 30 &  tempp{:,2} < 50),:);
   end
   
       tempp=  [Task{i,1}.peaks];
    if length(tempp) ==0
        gamma_Task{i,:,:}=[];
    else
    tempp=struct2table(tempp);
    gamma_Task{i,:,:}=tempp(find(tempp{:,2} >= 30 &  tempp{:,2} < 50),:);
    end
        
end


%% wrangle peak alpha into channel matrix 

count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(alpha_Task{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([alpha_Rest{i,1}.channel], RowNames(j)));

        end
    if ~ isnan(chanind)
        
        amplitude_alpha_Rest(count,j)=[alpha_Rest{i,1}(chanind,:).amplitude(1)];
        centralfreq_alpha_Rest(count,j)=[alpha_Rest{i,1}(chanind,:).center_frequency(1)];
        std_alpha_Rest(count,j)=[alpha_Rest{i,1}(chanind,:).std_dev(1)];
        alpha_peak_count_rest(count,j)=length(chanind);

    else
       alpha_peak_count_rest(count,j)=0;
       amplitude_alpha_Rest(count,j)=NaN;
       centralfreq_alpha_Rest(count,j)=NaN;
       std_alpha_Rest(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(alpha_Task{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([alpha_Task{i,1}.channel], RowNames(j)));
        end
    if ~ isnan(chanind)
        
        amplitude_alpha_Task(count,j)=[alpha_Task{i,1}(chanind,:).amplitude(1)];
        centralfreq_alpha_Task(count,j)=[alpha_Task{i,1}(chanind,:).center_frequency(1)];
        std_alpha_Task(count,j)=[alpha_Task{i,1}(chanind,:).std_dev(1)];
        alpha_peak_count_task(count,j)=length(chanind);

    else
       alpha_peak_count_task(count,j)=0;
       amplitude_alpha_Task(count,j)=NaN;
       centralfreq_alpha_Task(count,j)=NaN;
       std_alpha_Task(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(beta_Rest{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([beta_Rest{i,1}.channel], RowNames(j)));
        end
    
    if ~ isnan(chanind)
        
        amplitude_beta_Rest(count,j)=[beta_Rest{i,1}(chanind,:).amplitude(1)];
        centralfreq_beta_Rest(count,j)=[beta_Rest{i,1}(chanind,:).center_frequency(1)];
        std_beta_Rest(count,j)=[beta_Rest{i,1}(chanind,:).std_dev(1)];
        beta_peak_count_rest(count,j)=length(chanind);

    else
       beta_peak_count_rest(count,j)=0;
       amplitude_beta_Rest(count,j)=NaN;
       centralfreq_beta_Rest(count,j)=NaN;
       std_beta_Rest(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(beta_Task{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([beta_Task{i,1}.channel], RowNames(j)));
        end
    
    if ~ isnan(chanind)
        
        amplitude_beta_Task(count,j)=[beta_Task{i,1}(chanind,:).amplitude(1)];
        centralfreq_beta_Task(count,j)=[beta_Task{i,1}(chanind,:).center_frequency(1)];
        std_beta_Task(count,j)=[beta_Task{i,1}(chanind,:).std_dev(1)];
        beta_peak_count_task(count,j)=length(chanind);

    else
       beta_peak_count_task(count,j)=0;
       amplitude_beta_Task(count,j)=NaN;
       centralfreq_beta_Task(count,j)=NaN;
       std_beta_Task(count,j)=NaN;
        
    end
    end
    count=count+1;
end



% delta
count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(delta_Rest{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([delta_Rest{i,1}.channel], RowNames(j)));
        end
    
    if ~ isnan(chanind)
        
        amplitude_delta_Rest(count,j)=[delta_Rest{i,1}(chanind,:).amplitude(1)];
        centralfreq_delta_Rest(count,j)=[delta_Rest{i,1}(chanind,:).center_frequency(1)];
        std_delta_Rest(count,j)=[delta_Rest{i,1}(chanind,:).std_dev(1)];
        delta_peak_count_rest(count,j)=length(chanind);

    else
       delta_peak_count_rest(count,j)=0;
       amplitude_delta_Rest(count,j)=NaN;
       centralfreq_delta_Rest(count,j)=NaN;
       std_delta_Rest(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(delta_Task{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([delta_Task{i,1}.channel], RowNames(j)));
        end
    
    if ~ isnan(chanind)
        
        amplitude_delta_Task(count,j)=[delta_Task{i,1}(chanind,:).amplitude(1)];
        centralfreq_delta_Task(count,j)=[delta_Task{i,1}(chanind,:).center_frequency(1)];
        std_delta_Task(count,j)=[delta_Task{i,1}(chanind,:).std_dev(1)];
        delta_peak_count_task(count,j)=length(chanind);

    else
       delta_peak_count_task(count,j)=0;
       amplitude_delta_Task(count,j)=NaN;
       centralfreq_delta_Task(count,j)=NaN;
       std_delta_Task(count,j)=NaN;
        
    end
    end
    count=count+1;
end

% theta
count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(theta_Rest{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([theta_Rest{i,1}.channel], RowNames(j)));
        end
    
    if ~ isnan(chanind)
        
        amplitude_theta_Rest(count,j)=[theta_Rest{i,1}(chanind,:).amplitude(1)];
        centralfreq_theta_Rest(count,j)=[theta_Rest{i,1}(chanind,:).center_frequency(1)];
        std_theta_Rest(count,j)=[theta_Rest{i,1}(chanind,:).std_dev(1)];
        theta_peak_count_rest(count,j)=length(chanind);

    else
       theta_peak_count_rest(count,j)=0;
       amplitude_theta_Rest(count,j)=NaN;
       centralfreq_theta_Rest(count,j)=NaN;
       std_theta_Rest(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(theta_Task{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([theta_Task{i,1}.channel], RowNames(j)));
        end
    
    if ~ isnan(chanind)
        
        amplitude_theta_Task(count,j)=[theta_Task{i,1}(chanind,:).amplitude(1)];
        centralfreq_theta_Task(count,j)=[theta_Task{i,1}(chanind,:).center_frequency(1)];
        std_theta_Task(count,j)=[theta_Task{i,1}(chanind,:).std_dev(1)];
        theta_peak_count_task(count,j)=length(chanind);

    else
       theta_peak_count_task(count,j)=0;
       amplitude_theta_Task(count,j)=NaN;
       centralfreq_theta_Task(count,j)=NaN;
       std_theta_Task(count,j)=NaN;
        
    end
    end
    count=count+1;
end



% gamma
count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(gamma_Rest{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([gamma_Rest{i,1}.channel], RowNames(j)));
        end
    
    if ~ isnan(chanind)
        
        amplitude_gamma_Rest(count,j)=[gamma_Rest{i,1}(chanind,:).amplitude(1)];
        centralfreq_gamma_Rest(count,j)=[gamma_Rest{i,1}(chanind,:).center_frequency(1)];
        std_gamma_Rest(count,j)=[gamma_Rest{i,1}(chanind,:).std_dev(1)];
        gamma_peak_count_rest(count,j)=length(chanind);

    else
       gamma_peak_count_rest(count,j)=0;
       amplitude_gamma_Rest(count,j)=NaN;
       centralfreq_gamma_Rest(count,j)=NaN;
       std_gamma_Rest(count,j)=NaN;
        
    end
    end
    count=count+1;
end


count=1;
for i=1:numpar
    for j=1:nchan
        
        if height(gamma_Task{i,1}) ==0
            chanind=NaN;
        else
             chanind=find(strcmp([gamma_Task{i,1}.channel], RowNames(j)));
        end
    
    if ~ isnan(chanind)
        
        amplitude_gamma_Task(count,j)=[gamma_Task{i,1}(chanind,:).amplitude(1)];
        centralfreq_gamma_Task(count,j)=[gamma_Task{i,1}(chanind,:).center_frequency(1)];
        std_gamma_Task(count,j)=[gamma_Task{i,1}(chanind,:).std_dev(1)];
        gamma_peak_count_task(count,j)=length(chanind);

    else
       gamma_peak_count_task(count,j)=0;
       amplitude_gamma_Task(count,j)=NaN;
       centralfreq_gamma_Task(count,j)=NaN;
       std_gamma_Task(count,j)=NaN;
        
    end
    end
    count=count+1;
end


spectral_peak_count_rest = [sum(delta_peak_count_rest,1); sum(theta_peak_count_rest,1);...
    sum(alpha_peak_count_rest,1); sum(beta_peak_count_rest,1);sum(gamma_peak_count_rest,1)];


spectral_peak_count_task = [sum(delta_peak_count_task,1); sum(theta_peak_count_task,1);...
    sum(alpha_peak_count_task,1); sum(beta_peak_count_task,1);sum(gamma_peak_count_task,1)];




%% Collect fits of aperiodic and periodic 


for i=1:numpar
    
    temp11=reshape([Rest{i,1}.data.ap_fit], [size([Rest{i,1}.data.ap_fit],2)/68,68]);
    temp12=reshape([Task{i,1}.data.ap_fit], [size([Task{i,1}.data.ap_fit],2)/68,68]);

    fitaper_Rest(i,:,:)=  temp11(1:78,:);
    fitaper_Task(i,:,:)=  temp12(1:78,:);
    
    temp11=reshape([Rest{i,1}.data.peak_fit], [size([Rest{i,1}.data.peak_fit],2)/68,68]);
    temp12=reshape([Task{i,1}.data.peak_fit], [size([Task{i,1}.data.peak_fit],2)/68,68]);
 
    fitper_Rest(i,:,:)=  temp11(1:78,:);
    fitper_Task(i,:,:)=  temp12(1:78,:);
 
end


Rest_psd=squeeze(log(mean(mean((fitper_Rest),3),1)));
Task_psd=squeeze(log(mean(mean((fitper_Task),3),1)));
figure
clear g

g(1,1)=gramm('x',repmat(Freqs(3:80),[2,1]),'y', [Rest_psd',Task_psd']', 'color',{'Rest', 'Task'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Periodic --Baseline');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();

Rest_psd=squeeze(log(mean(mean((fitaper_Rest),3),1)));
Task_psd=squeeze(log(mean(mean((fitaper_Task),3),1)));
figure
clear g

g(1,1)=gramm('x',repmat(Freqs(3:80),[2,1]),'y', [Rest_psd',Task_psd']', 'color',{'Rest', 'Task'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Aperiodic --Baseline');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();


gauss_psd=squeeze(mean(mean(log(fitaper_Rest),3),1));
cau_psd=squeeze(mean(mean(log(PSD_rest(:,:,3:80)),2),1));
figure
clear g

g(1,1)=gramm('x',repmat(log(Freqs(3:80)),[2,1]),'y', [gauss_psd',cau_psd]', 'color',{'FOOOF', 'PSD'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Fitted Spectrum --Rest');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();

saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_fitted_specra_rest.pdf')



gauss_psd=squeeze(mean(mean(log(fitaper_Task),3),1));
cau_psd=squeeze(mean(mean(log(PSD_task(:,:,3:80)),2),1));
figure
clear g

g(1,1)=gramm('x',repmat(log(Freqs(3:80)),[2,1]),'y', [gauss_psd',cau_psd]', 'color',{'FOOOF', 'PSD'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Fitted Spectrum --Task');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();

saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_fitted_specra_task.pdf')

%%
figure
histogram(offest_Rest)
hold on;
histogram(offest_Task)

%% Output

cd('/Users/jason/Documents/CAMCAN_outputs/CAMCAN_Analysis/')
writetable(cell2table(Ids), 'subject_codes.csv')
writecell(RowNames, 'ROIs.csv')

writematrix(expon_Rest, 'aperiodic_exponent_Rest.csv')
writematrix(expon_Task, 'aperiodic_exponent_Task.csv')
writematrix(offest_Rest, 'aperiodic_offset_Rest.csv')
writematrix(offest_Task, 'aperiodic_offset_Task.csv')

writematrix(mse_Rest, 'fit_mse_Rest.csv')
writematrix(mse_Task, 'fit_mse_Task.csv')
writematrix(rsq_Rest, 'fit_rsq_Rest.csv')
writematrix(rsq_Task, 'fit_rsq_Task.csv')

writematrix(amplitude_alpha_Rest, 'alpha_amplitude_Rest.csv')
writematrix(amplitude_alpha_Task, 'alpha_amplitude_Task.csv')
writematrix(amplitude_beta_Rest, 'beta_amplitude_Rest.csv')
writematrix(amplitude_beta_Task, 'beta_amplitude_Task.csv')

writematrix(centralfreq_alpha_Rest, 'alpha_centerfreq_Rest.csv')
writematrix(centralfreq_alpha_Task, 'alpha_centerfreq_Task.csv')
writematrix(centralfreq_beta_Rest, 'beta_centerfreq_Rest.csv')
writematrix(centralfreq_beta_Task, 'beta_centerfreq_Task.csv')

writematrix(std_alpha_Rest, 'alpha_std_Rest.csv')
writematrix(std_alpha_Task, 'alpha_std_Task.csv')
writematrix(std_beta_Rest, 'beta_std_Rest.csv')
writematrix(std_beta_Task, 'beta_std_Task.csv')

corr_spec_rest=log(permute(PSD_rest(:,:,3:80), [1,3,2]))-log(fitaper_Rest);
corr_spec_task=log(permute(PSD_task(:,:,3:80), [1,3,2]))-log(fitaper_Task);

figure
clear g

g(1,1)=gramm('x',Freqs(3:80),'y', squeeze(mean(corr_spec_rest,1))', 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Aperiodic corrected spectra Rest');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_aperiodic_corrected_rest.pdf')


figure
clear g

g(1,1)=gramm('x',Freqs(3:80),'y', squeeze(mean(corr_spec_task,1))', 'color', categorical(1:68));
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Aperiodic corrected spectra Task');
g(1,1).no_legend();
g(1,1).set_color_options('map',colourmap);
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();
saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_aperiodic_corrected_task.pdf')


gauss_psd=squeeze(mean(mean(corr_spec_rest,3),1));
cau_psd=squeeze(mean(mean(corr_spec_task,3),1));;
figure
clear g

g(1,1)=gramm('x',repmat(log(Freqs(3:80)),[2,1]),'y', [gauss_psd;cau_psd], 'color',{'Rest', 'Task'});
%smooth plot 
g(1,1).geom_line();
g(1,1).set_title('Aperiodic corrected spectra');
g(1,1).set_color_options('map','brewer2');
%These functions can be called on arrays of gramm objects
g.set_names('x','Frequency (Hz)','y','Log Power');
%g.axe_property('YLim',[1e-15 1.5e-9]);
%g.axe_property('XLim',[0 50]);
g.draw();

saveas(gcf,'/Users/jason/Documents/CAMCAN_outputs/figures/MATLAB_aperiodic_corrected_spectra_mean.pdf')


writematrix(reshape(corr_spec_rest, [606,78*68]), 'corrspec_Rest.csv')
writematrix(reshape(corr_spec_task, [606,78*68]), 'corrspec_Task.csv')


figure
plot(squeeze(mean(corr_spec_rest,1)))

writematrix(reshape(permute(PSD_rest(:,:,1:300), [1,3,2]), [606,68*300]), 'PSD_Rest.csv')
writematrix(reshape(permute(PSD_task(:,:,1:300), [1,3,2]), [606,68*300]), 'PSD_Task.csv')

writematrix(reshape(fitaper_Rest, [606,68*78]), 'aperioidc_Rest.csv')
writematrix(reshape(fitaper_Task, [606,68*78]), 'aperiodic_Task.csv')

writematrix(reshape(permute(fitper_Rest, [1,3,2]), [606,68*78]), 'perioidc_Rest.csv')
writematrix(reshape(permute(fitper_Task, [1,3,2]), [606,68*78]), 'periodic_Task.csv')


writematrix(spectral_peak_count_rest', 'spectral_peak_count_Rest.csv')
writematrix(spectral_peak_count_task', 'spectral_peak_count_Task.csv')


writematrix(delta_peak_count_rest, 'spectral_peak_count_rest_delta.csv')
writematrix(theta_peak_count_rest, 'spectral_peak_count_rest_theta.csv')
writematrix(alpha_peak_count_rest, 'spectral_peak_count_rest_alpha.csv')
writematrix(beta_peak_count_rest, 'spectral_peak_count_rest_beta.csv')
writematrix(gamma_peak_count_rest, 'spectral_peak_count_rest_gamma.csv')

%%

cd('/Users/jason/Documents/CAMCAN_outputs/fingerprint_data/')

for i=1:numpar
    temp=squeeze(PSD_rest(i,:,:));
    temp2=squeeze(PSD_task(i,:,:));
   
    csvwrite( strcat('spectra/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('spectra/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
    temp=squeeze(fitaper_Rest(i,:,:))';
    temp2=squeeze(fitaper_Task(i,:,:))';
    
    csvwrite( strcat('FOOOF_aperiodic/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('FOOOF_aperiodic/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
    temp=squeeze(fitper_Rest(i,:,:))';
    temp2=squeeze(fitper_Task(i,:,:))';
    
    csvwrite( strcat('FOOOF_periodic/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('FOOOF_periodic/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
    
    temp=squeeze(fmse_Rest(i,:,:))';
    temp2=squeeze(fmse_Task(i,:,:))';
    
    csvwrite( strcat('FOOOF_residuals/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('FOOOF_residuals/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
    
    temp=squeeze(fpsd_Rest(i,:,:))';
    temp2=squeeze(fpsd_Task(i,:,:))';
    
    csvwrite( strcat('FOOOF_residuals/sub_', num2str(i, '%04.f'), '_training.csv'), temp);
    csvwrite(strcat('FOOOF_residuals/sub_', num2str(i,'%04.f'), '_validation.csv'), temp2);
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%               Process just baseline
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

colourmap=repmat([90:2.2:255],[3,1])'./255;
nchan=68;
numpar=606;
% Wrangle data
files_RS=dir('/Users/jason/Documents/CAMCAN_outputs/CAMCAN_just_baseline_psd/*.mat');

count1=1;
count2=1;
for i=1:numpar*2
    
    if mod(i,2)~=0
        load(strcat(files_RS(i,:).folder,'/',files_RS(i,:).name));
        PSD_rest(count1,:,:)=squeeze(TF);
        count1=count1+1;
    
    else
    
        load(strcat(files_RS(i,:).folder,'/',files_RS(i,:).name));
        PSD_task(count2,:,:)=squeeze(TF);
        count2=count2+1;

    end
    
    
end 

cd('/Users/jason/Documents/CAMCAN_outputs/CAMCAN_Analysis/')


writematrix(reshape(permute(PSD_rest(:,:,1:300), [1,3,2]), [606,68*300]), 'PSD_Rest1.csv')
writematrix(reshape(permute(PSD_task(:,:,1:300), [1,3,2]), [606,68*300]), 'PSD_Rest2.csv')
