%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to calculate the mean and standard errors of the simulation
%combinations for the Bayesian nonparanormal method
% Author: Jami Jackson Mulgrave
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%AR1 model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p=100 n=500

load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_1to10_final.mat');


%Vertically concatenate the row vectors
SE_n500_p100_AR1_vector = zeros([reps,1]); 
SP_n500_p100_AR1_vector = zeros([reps,1]); 
MCC_n500_p100_AR1_vector = zeros([reps,1]);
total_time_n500_p100_AR1_vector = zeros([reps,1]);

SE_n500_p100_AR1_EBIC_vector = zeros([reps,1]); 
SP_n500_p100_AR1_EBIC_vector = zeros([reps,1]); 
MCC_n500_p100_EBIC_vector = zeros([reps,1]);



for iters = 1:10
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;

end

%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis

load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_11to20_final.mat');


for iters = 11:20
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;

end

%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis


load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_21to30_final.mat');


for iters = 21:30
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;

end

%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis

load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_31to40_final.mat');

for iters = 31:40
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;
end


%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis


load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_41to50_final.mat');

for iters = 41:50
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;
end



%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis


load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_51to60_final.mat');

for iters = 51:60
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;
end


%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis


load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_61to70_final.mat');

for iters = 61:70
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;
end


%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis


load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_71to80_final.mat');

for iters = 71:80
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;
end


%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis


load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_81to90_final.mat');

for iters = 81:90
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;
end


%Clear these matrices because they are named the same in the next datafile
clear total_time_n500_p100_AR1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;
clear SP_matrix_eBIC_finalanalysis;
clear SE_matrix_eBIC_finalanalysis;
clear MCC_matrix_eBIC_finalanalysis


load('BayesNonpar_Sim_n500_p100_AR1_SS_1to100_91to100_final.mat');

for iters = 91:100
    total_time_n500_p100_AR1_vector(iters) = total_time_n500_p100_AR1_finalanalysis(iters);
    SP_n500_p100_AR1_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n500_p100_AR1_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n500_p100_AR1_vector(iters) =   MCC_matrix_finalanalysis(iters);
    
   SP_n500_p100_AR1_EBIC_vector(iters) =  SP_matrix_eBIC_finalanalysis(iters);
   SE_n500_p100_AR1_EBIC_vector(iters) = SE_matrix_eBIC_finalanalysis(iters);
   MCC_n500_p100_EBIC_vector(iters) =  MCC_matrix_eBIC_finalanalysis(iters) ;
end


%find the mean and std errors

total_time_n500_p100_AR1_mean_rounded = round2(mean(total_time_n500_p100_AR1_vector),.01);
total_time_n500_p100_AR1_SE_rounded = round2(std(total_time_n500_p100_AR1_vector)/...
    sqrt(length(total_time_n500_p100_AR1_vector)),.01);
 
    
     SP_matrix_eBIC_mean_rounded = round2(mean(SP_n500_p100_AR1_EBIC_vector),.01);
     SP_matrix_eBIC_se_rounded = round2(std(SP_n500_p100_AR1_EBIC_vector)/...
         sqrt(length(SP_n500_p100_AR1_EBIC_vector)), .01);
     
    SE_matrix_eBIC_mean_rounded = round2(mean(SE_n500_p100_AR1_EBIC_vector), .01);
    SE_matrix_eBIC_se_rounded = round2(std(SE_n500_p100_AR1_EBIC_vector)/...
        sqrt(length(SE_n500_p100_AR1_EBIC_vector)), .01);
    
    MCC_matrix_eBIC_mean_rounded = round2(nanmean(MCC_n500_p100_EBIC_vector), .01);
    MCC_matrix_eBIC_se_rounded = round2(nanstd(MCC_n500_p100_EBIC_vector)/...
        sqrt(length(MCC_n500_p100_EBIC_vector(~isnan(MCC_n500_p100_EBIC_vector)))),.01);
   
%find the mean and SE for  SE, SP, MCC

% I have to use round2 and use .01 for the number of significant digits

SE_n500_p100_AR1_mean_rounded = round2(mean(SE_n500_p100_AR1_vector),.01);
SE_n500_p100_AR1_SE_rounded = round2(std(SE_n500_p100_AR1_vector)/...
    sqrt(length(SE_n500_p100_AR1_vector)),.01);

SP_n500_p100_AR1_mean_rounded = round2(mean(SP_n500_p100_AR1_vector),.01);
SP_n500_p100_AR1_SE_rounded = round2(std(SP_n500_p100_AR1_vector)/...
    sqrt(length(SP_n500_p100_AR1_vector)),.01);

%have to exclude NaNs when reporting MCC.
MCC_n500_p100_AR1_mean_rounded = round2(nanmean(MCC_n500_p100_AR1_vector),.01);
MCC_n500_p100_AR1_SE_rounded = round2(nanstd(MCC_n500_p100_AR1_vector)/...
    sqrt(length(MCC_n500_p100_AR1_vector(~isnan(MCC_n500_p100_AR1_vector)))),.01);


save('Bsplinespaper_Calcs_p100_n500_AR1_SpikeSlab_1to100.mat',...
    'total_time_n500_p100_AR1_mean_rounded',... 
'total_time_n500_p100_AR1_SE_rounded',...
   'SE_n500_p100_AR1_mean_rounded',...
'SE_n500_p100_AR1_SE_rounded', 'SP_n500_p100_AR1_mean_rounded', 'SP_n500_p100_AR1_SE_rounded',...
'MCC_n500_p100_AR1_mean_rounded',...
 'MCC_n500_p100_AR1_SE_rounded',...
 'SP_matrix_eBIC_mean_rounded', 'SP_matrix_eBIC_se_rounded', 'SE_matrix_eBIC_mean_rounded',...
 'SE_matrix_eBIC_se_rounded',...
 'MCC_matrix_eBIC_mean_rounded', 'MCC_matrix_eBIC_se_rounded');

