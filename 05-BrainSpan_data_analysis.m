% This script organizes and parcellates data from BrainSpan

% BrainSpan data can be downloaded in its original form from
% https://www.brainspan.org/static/download.html

%% load
clear all
close all

cd '~/Documents/HCP_twin_projec/abagen_analysis/'
load('rowids.mat') % node by gene expression matrix
load('gene_expression_abagen.mat')           % relevant node indices
load('gene_names.mat')          
load('result_bandlimited.mat')
load('SPINStwirls.mat')           % spin test indices

label=genenames.genenames; % genes in AHBA

spins= permutedindexesofschaeferatlasSPINsTwirl+1;


expressiongenes= expressiongenes(rowids,:); % fix row ids

cd('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/hansen_genescognition-master');
load('./BrainSpan/mapping.mat')  % mapping from 34 node parcellation to the 16 unique cortical regions included in BrainSpan

% AHBA harmonized files have been organized to include only genes included
% in AHBA
brainspan = readtable('./BrainSpan/gene_expression_AHBA_harmonized.csv'); % gene by sample matrix of gene expression
gene_info = readtable('./BrainSpan/gene_metadata_AHBA_harmonized.csv');   % metadata on genes
sample_info = readtable('./BrainSpan/samples_metadata.csv');              % metadata on tissue samples
brainspan = table2array(brainspan(2:end,2:end));              % remove column of gene names and row of sample IDs

%% remove non-cortical samples
% remove samples labeled 'not cortex' and amygdala samples

notcortex_idx = [];
sensorifugal = table2cell(sample_info(:,15)); % this column labels noncortical regions as 'Not_Cortex' 
sensorifugal = string(sensorifugal);

% find noncortical defined by sensorifugal
for k = 1:length(sensorifugal)
    if strcmp(sensorifugal(k),'Not_Cortex')
        notcortex_idx = [notcortex_idx; k];
    end
end

% find amygdala indices
amygdala = find(contains(table2cell(sample_info(:,8)),'amygdaloid complex'));
notcortex_idx = [notcortex_idx; amygdala];
notcortex_idx = sort(notcortex_idx);

% get all relevant (cortical) indices
cortex_idx = setdiff([1:length(sensorifugal)],notcortex_idx); 

% remove noncortical indices
brainspan(:,notcortex_idx) = [];


%% find pos neg loadings

% get gene sets
% compute the loading of each gene as the correlation between the original
% data and the gene scores
ngenes = length(expressiongenes);
gload = zeros(ngenes,1);
for k = 1:ngenes
    gload(k) = corr(expressiongenes(:,k),result.vsc(:,1));
end

%gload2 = corr(expressiongenes,result.vsc(:,1)); Does the same thing 

ipos = find(gload > 0); % index of genes with positive loading
ineg = find(gload < 0); % index of genes with negative loading
gload_pos = gload(gload > 0); % loading of genes with positive loading
gload_neg = gload(gload < 0); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos); % sorted
[~,Ineg] = sort(gload_neg); % sorted 

threshold = 0.5; % top 50% of pos/neg genes constitute each gene set

gpos_idx = Ipos(end-floor(threshold*length(gload_pos)):end); % top 50% of genes with positive loading
gneg_idx = Ineg(1:floor(threshold*length(gload_neg)));       % top 50% of genes with negative loading


%% get index of pos genes

% genes included in original analyses are stable (differential stability >
% 0.1). For comparability, we only include (available) stable genes as defined on the
% our schaefer parcellation

pos_label = label(ipos(gpos_idx)); % names of positive loaded genes
gene_label = table2cell(gene_info(:,4));     % names of genes in BrainSpan

notstable_idx = [];
bspan_gidx033 = [];
for k = 1:length(gene_label)                            % for each gene in BrainSpan
    if ~ismember(gene_label(k), pos_label)            % if gene isn't stable
        notstable_idx = [notstable_idx; k];             % keep its index
    else                                                % if gene is stable
        i = find(ismember(pos_label,gene_label(k))); % get index of stable gene
        bspan_gidx033 = [bspan_gidx033; i];             % store it
    end
end

% remove nonstable genes
pos_brainspan= brainspan;
pos_brainspan(notstable_idx,:) = [];


neg_label = label(ineg(gneg_idx)); % names of positive loaded genes
gene_label = table2cell(gene_info(:,4));     % names of genes in BrainSpan

notstable_idx = [];
bspan_gidx033 = [];
for k = 1:length(gene_label)                            % for each gene in BrainSpan
    if ~ismember(gene_label(k), neg_label)            % if gene isn't stable
        notstable_idx = [notstable_idx; k];             % keep its index
    else                                                % if gene is stable
        i = find(ismember(neg_label,gene_label(k))); % get index of stable gene
        bspan_gidx033 = [bspan_gidx033; i];             % store it
    end
end

% remove nonstable genes
neg_brainspan= brainspan;
neg_brainspan(notstable_idx,:) = [];


%% organize data by life stage
% samples are organized into 5 life stages: fetal, infant, child,
% adolescent, and adult

% get indices
fetal_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'fetal'));
infant_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'infant'));
child_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'child'));
adolescent_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adolescent'));
adult_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adult'));

% make gene expression matrices for each life stage
fetal = neg_brainspan(:,fetal_idx);
infant = neg_brainspan(:,infant_idx);
child = neg_brainspan(:,child_idx);
adolescent = neg_brainspan(:,adolescent_idx);
adult = neg_brainspan(:,adult_idx);

% organize matrices and indices
M_idx = {fetal_idx, infant_idx, child_idx, adolescent_idx, adult_idx};
M_NEG = {fetal, infant, child, adolescent, adult};

fetal = pos_brainspan(:,fetal_idx);
infant = pos_brainspan(:,infant_idx);
child = pos_brainspan(:,child_idx);
adolescent = pos_brainspan(:,adolescent_idx);
adult = pos_brainspan(:,adult_idx);

% organize matrices and indices
M_idx = {fetal_idx, infant_idx, child_idx, adolescent_idx, adult_idx};
M_POS = {fetal, infant, child, adolescent, adult};

% parcellate brain regions

sample_regions = table2cell(sample_info(cortex_idx,8)); % get region name of each sample
regions = unique(string(sample_regions));               % get unique regions - this is the maximum number of brain regions in each gene expression matrix

% make value-based region mapping instead of string-based (for
% simplicity)
region_mapping = zeros(length(sample_regions),1);

for k = 1:length(regions)                               % for each unique region
    i = find(contains(sample_regions,regions(k)));      % find all samples from that region
    region_mapping(i) = k;                              % make index-based map of regions
end

% average expression of each gene from identical regions
% done separately for each life stage, where the number of unique regions
% with gene expression estimates varies across life stage
regionsIncluded = cell(5,1);
for k = 1:length(M_NEG)                                          % for each life stage
    mat = M_NEG{k};                                              % get original gene expression matrix (genes x samples)
    r = region_mapping(M_idx{k});                            % get region mapping of this life stage
    regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
    i = 1;
    for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
        i = i+1;
    end
    M_NEG{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
end

regionsIncluded = cell(5,1);
for k = 1:length(M_POS)                                          % for each life stage
    mat = M_POS{k};                                              % get original gene expression matrix (genes x samples)
    r = region_mapping(M_idx{k});                            % get region mapping of this life stage
    regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
    i = 1;
    for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
        i = i+1;
    end
    M_POS{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
end


% remove the regions only present in feotus
regionsindfeotus= [1:4, 6, 8:11, 13, 14, 16];
M_POS{1}= M_POS{1}(regionsindfeotus,:);
M_NEG{1}= M_NEG{1}(regionsindfeotus,:);

negative_topo= zeros([12,5]);
positive_topo= zeros([12,5]);

for k = 1:length(M_POS)       
    
    temp= M_POS{k};
    positive_topo(:,k) = mean(temp,2, 'omitnan');
    
    temp= M_NEG{k};
    negative_topo(:,k) = mean(temp,2, 'omitnan');

end

figure
r=corr(positive_topo)
isupper = logical(triu(ones(size(r)),1));
r(isupper) = NaN;
% Plot results
h = heatmap(r,'MissingDataColor','w');
labels = ["fetal","infant","child","adolescent","adult"];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels; 
colormap (flip(jet));

figure
r=corr(negative_topo)
isupper = logical(triu(ones(size(r)),1));
r(isupper) = NaN;
% Plot results
h = heatmap(r,'MissingDataColor','w');
labels = ["fetal","infant","child","adolescent","adult"];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels; 
colormap (flip(jet));


reorderind=[5,12,1,2,4,3,8,6,9,10,7,11];
regionsfinal=regions(regionsindfeotus)
figure
scatter(repelem(1:12,5)',reshape(negative_topo(reorderind,:)',1,[]), 100,'filled')
hold on;
plot(negative_topo(reorderind,:), 'linewidth',2)
legend
labels = regionsfinal(reorderind);
xticks(1:12)
xticklabels(labels)


%% change definition of child

% change child sample to exclude 2 and 3 yrs old to match MEG
sample_info{strcmp(sample_info{:,4}, '2 yrs'),9} = repelem({'toddler'}, length(find(strcmp(sample_info{:,4}, '2 yrs'))))';
sample_info{strcmp(sample_info{:,4}, '3 yrs'),9} = repelem({'toddler'}, length(find(strcmp(sample_info{:,4}, '3 yrs'))))';


% get indices
fetal_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'fetal'));
infant_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'infant'));
toddler_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'toddler'));
child_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'child'));
adolescent_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adolescent'));
adult_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adult'));

% make gene expression matrices for each life stage
fetal = neg_brainspan(:,fetal_idx);
infant = neg_brainspan(:,infant_idx);
toddler = neg_brainspan(:,toddler_idx);
child = neg_brainspan(:,child_idx);
adolescent = neg_brainspan(:,adolescent_idx);
adult = neg_brainspan(:,adult_idx);

% organize matrices and indices
M_idx = {fetal_idx, infant_idx, toddler_idx,child_idx, adolescent_idx, adult_idx};
M_NEG = {fetal, infant, toddler, child, adolescent, adult};

fetal = pos_brainspan(:,fetal_idx);
infant = pos_brainspan(:,infant_idx);
toddler = pos_brainspan(:,toddler_idx);
child = pos_brainspan(:,child_idx);
adolescent = pos_brainspan(:,adolescent_idx);
adult = pos_brainspan(:,adult_idx);

% organize matrices and indices
M_idx = {fetal_idx, infant_idx, toddler_idx,child_idx, adolescent_idx, adult_idx};
M_POS = {fetal, infant, toddler, child, adolescent, adult};

% parcellate brain regions

sample_regions = table2cell(sample_info(cortex_idx,8)); % get region name of each sample
regions = unique(string(sample_regions));               % get unique regions - this is the maximum number of brain regions in each gene expression matrix

% make value-based region mapping instead of string-based (for
% simplicity)
region_mapping = zeros(length(sample_regions),1);

for k = 1:length(regions)                               % for each unique region
    i = find(contains(sample_regions,regions(k)));      % find all samples from that region
    region_mapping(i) = k;                              % make index-based map of regions
end

% average expression of each gene from identical regions
% done separately for each life stage, where the number of unique regions
% with gene expression estimates varies across life stage
regionsIncluded = cell(5,1);
for k = 1:length(M_NEG)                                          % for each life stage
    mat = M_NEG{k};                                              % get original gene expression matrix (genes x samples)
    r = region_mapping(M_idx{k});                            % get region mapping of this life stage
    regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
    i = 1;
    for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
        i = i+1;
    end
    M_NEG{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
end

regionsIncluded = cell(5,1);
for k = 1:length(M_POS)                                          % for each life stage
    mat = M_POS{k};                                              % get original gene expression matrix (genes x samples)
    r = region_mapping(M_idx{k});                            % get region mapping of this life stage
    regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
    i = 1;
    for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
        i = i+1;
    end
    M_POS{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
end


% remove the regions only present in feotus
regionsindfeotus= [1:4, 6, 8:11, 13, 14, 16];
M_POS{1}= M_POS{1}(regionsindfeotus,:);
M_NEG{1}= M_NEG{1}(regionsindfeotus,:);

negative_topo= zeros([12,6]);
positive_topo= zeros([12,6]);

for k = 1:length(M_POS)       
    
    temp= M_POS{k};
    positive_topo(:,k) = mean(temp,2, 'omitnan');
    
    temp= M_NEG{k};
    negative_topo(:,k) = mean(temp,2, 'omitnan');

end

figure
r=corr(positive_topo)
isupper = logical(triu(ones(size(r)),1));
r(isupper) = NaN;
% Plot results
h = heatmap(r,'MissingDataColor','w');
labels = ["fetal","infant",'toddler',"child","adolescent","adult"];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels; 
colormap (flip(jet));

figure
r=corr(negative_topo)
isupper = logical(triu(ones(size(r)),1));
r(isupper) = NaN;
% Plot results
h = heatmap(r,'MissingDataColor','w');
labels = ["fetal","infant",'toddler',"child","adolescent","adult"];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels; 
colormap (flip(jet));



reorderind=[5,12,1,2,4,3,8,6,9,10,7,11];
regionsfinal=regions(regionsindfeotus)
figure
scatter(negative_topo(reorderind,:)')
hold on;
plot(negative_topo(reorderind,:), 'linewidth',2)
legend
labels = regionsfinal(reorderind);
xticks(1:12)
xticklabels(labels)

