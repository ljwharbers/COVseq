%This script was used to calculate the cumulative costs of making SARS-CoV-2 WGS libraries using
%COVseq vs. three commercial kits, as summarized in Figure 4 and described
%in details in the Supplementary Notes.
%
%Written by Nicola Crosetto in January 2021 using MATLAB version R2018b.

close all
clear all

%--------------------------------------------------------------------------
%GENERAL CONSIDERATIONS
%--------------------------------------------------------------------------
%In order to run, this script must be placed in the same folder as the
%Excel files provided in the same Github repository as this file.
%These files contain information about the price (in USD) and
%number of reactions of each reagent used.
%
%The script performs a cumulative cost analysis for a given number of
%samples processed, updating the cumulative costs whenever a reagent is
%finished and must be purchased again.
%
%NOTE: the script only calculates the cumulative costs for making
%libraries and does not include the costs for sequencing the samples
%(these depend on the type of Illumina platform used, on the
%sequencing mode, on the number of samples, and on the target sequencing
%depth per sample).

%--------------------------------------------------------------------------
%SET PARAMETERS
%--------------------------------------------------------------------------
%Specify the path where this file and the associated Excel files are located on your computer
pathname='';
%Set cost analysis parameters
tot_samp=10000; %Total # of samples
multiplex_num=96; %# of barcodes in COVseq or indexes in commercial kits
%NOTE:the analysis shown in Fig. 4 was done assuming using 96 COVseq
%barcodes and 96 indexes for the commercial kits

%--------------------------------------------------------------------------
%CALCULATE CUMULATIVE COSTS
%--------------------------------------------------------------------------
methods={'covseq-cdc','covseq-artic','cleanplex','nebnext','nextera'};
for k=1:numel(methods) %Loop through each method
    %Load data
    filename{k,1}=[pathname,'/',methods{k},'_reagent_costs_',num2str(multiplex_num),'.xlsx'];
    reagents=xlsread(filename{k,1});
    condition{k,1}=[methods{k},'-',num2str(multiplex_num),'-',num2str(tot_samp)];
    cost_per_sample(1,1)=sum(reagents(1:end-1,1));
    cum_cost(1,k)=cost_per_sample;
    if isempty(strfind(methods{k},'covseq'))==0
        lib_num=0; %Counter of libraries
        for i=2:tot_samp
            sample_cost_idx=find(reagents(:,3)==1);
            lib_cost_idx=find(reagents(:,4)==1);
            seq_cost_idx=find(reagents(:,5)==1);
            %Update costs for RT, multiplexed PCR and CUTseq
            for j=1:numel(sample_cost_idx) %Loop through each cost for sample
                if rem(i,reagents(sample_cost_idx(j),2))==0                    
                    cost_per_sample=cost_per_sample+reagents(sample_cost_idx(j),1);
                end
            end    
            %Update costs for library prep
            if rem(i,96)==0
                for j=1:numel(lib_cost_idx) %Loop through each cost for library
                    if rem(lib_num,reagents(lib_cost_idx(j),2))==0                    
                        cost_per_sample=cost_per_sample+reagents(lib_cost_idx(j),1);
                    end 
                end
                lib_num=lib_num+1;
            end
            cum_cost(i,k)=cost_per_sample;
        end
    else
        for i=2:tot_samp %1 sample = 1 library
            lib_cost_idx=find(reagents(:,3)==1);
            seq_cost_idx=find(reagents(:,4)==1);  
            %Update costs for library prep
            for j=1:numel(lib_cost_idx) %Loop through each cost for library
                if rem(i,reagents(lib_cost_idx(j),2))==0                    
                    cost_per_sample=cost_per_sample+reagents(lib_cost_idx(j),1); 
                end
            end
            cum_cost(i,k)=cost_per_sample;
        end
    end
end
mean_sample_cost=cum_cost(end,:)/tot_samp;

%--------------------------------------------------------------------------
%PLOT CUMULATIVE COST CURVES
%--------------------------------------------------------------------------
figure(1)
plot(1:tot_samp,cum_cost)
legend(condition);
figure(2)
bar(mean_sample_cost);
set(gca,'xticklabel',condition)
