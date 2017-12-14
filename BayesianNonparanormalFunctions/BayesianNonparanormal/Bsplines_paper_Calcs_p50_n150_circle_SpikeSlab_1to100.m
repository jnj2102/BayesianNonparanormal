%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to calculate the mean and standard errors of the simulation
%combinations for the Bayesian nonparanormal method
% Author: Jami Jackson Mulgrave
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%circle model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p=50 n=150

load('BayesNonpar_Sim_n150_p50_circle_SS_1to100_1to25_final.mat');


%Vertically concatenate the row vectors of the Bayes L1 loss 
SE_n150_p50_circle_vector = zeros([reps,1]); 
SP_n150_p50_circle_vector = zeros([reps,1]); 
MCC_n150_p50_circle_vector = zeros([reps,1]);
L1_n150_p50_circle_vector = zeros([reps,1]);
total_time_n150_p50_circle_vector = zeros([reps,1]);



for iters = 1:25
    total_time_n150_p50_circle_vector(iters) = total_time_n150_p50_circle_finalanalysis(iters);
    L1_n150_p50_circle_vector(iters) = L1_finalanalysis(iters);
    SP_n150_p50_circle_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n150_p50_circle_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n150_p50_circle_vector(iters) =   MCC_matrix_finalanalysis(iters);

end

%Clear these matrices because they are named the same in the next datafile
clear total_time_n150_p50_circle_finalanalysis;
clear L1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;

load('BayesNonpar_Sim_n150_p50_circle_SS_1to100_26to50_final.mat');


for iters = 26:50
    total_time_n150_p50_circle_vector(iters) = total_time_n150_p50_circle_finalanalysis(iters);
    L1_n150_p50_circle_vector(iters) = L1_finalanalysis(iters);
    SP_n150_p50_circle_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n150_p50_circle_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n150_p50_circle_vector(iters) =   MCC_matrix_finalanalysis(iters);

end

clear total_time_n150_p50_circle_finalanalysis;
clear L1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;


load('BayesNonpar_Sim_n150_p50_circle_SS_1to100_51to75_final.mat');


for iters = 51:75
    total_time_n150_p50_circle_vector(iters) = total_time_n150_p50_circle_finalanalysis(iters);
    L1_n150_p50_circle_vector(iters) = L1_finalanalysis(iters);
    SP_n150_p50_circle_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n150_p50_circle_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n150_p50_circle_vector(iters) =   MCC_matrix_finalanalysis(iters);

end

clear total_time_n150_p50_circle_finalanalysis;
clear L1_finalanalysis;
clear SP_matrix_finalanalysis;
clear SE_matrix_finalanalysis;
clear MCC_matrix_finalanalysis;

load('BayesNonpar_Sim_n150_p50_circle_SS_1to100_76to100_final.mat');

for iters = 76:100
    total_time_n150_p50_circle_vector(iters) = total_time_n150_p50_circle_finalanalysis(iters);
    L1_n150_p50_circle_vector(iters) = L1_finalanalysis(iters);
    SP_n150_p50_circle_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n150_p50_circle_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n150_p50_circle_vector(iters) =   MCC_matrix_finalanalysis(iters);

end


%find the mean time and L1

total_time_n150_p50_circle_mean_rounded = round2(mean(total_time_n150_p50_circle_vector),.01);
total_time_n150_p50_circle_SE_rounded = round2(std(total_time_n150_p50_circle_vector)/sqrt(length(total_time_n150_p50_circle_vector)),.01);


L1_n150_p50_circle_mean_rounded = round2(mean(L1_n150_p50_circle_vector),.01);
L1_n150_p50_circle_SE_rounded = round2(std(L1_n150_p50_circle_vector)/sqrt(length(L1_n150_p50_circle_vector)),.01);



%find the mean and SE for  SE, SP, MCC

%Since the cluster has an older version of Matlab, the round2 function is different.  I have to use round2 and use .01 for the number of significant digits

SE_n150_p50_circle_mean_rounded = round2(mean(SE_n150_p50_circle_vector),.01);
SE_n150_p50_circle_SE_rounded = round2(std(SE_n150_p50_circle_vector)/sqrt(length(SE_n150_p50_circle_vector)),.01);

SP_n150_p50_circle_mean_rounded = round2(mean(SP_n150_p50_circle_vector),.01);
SP_n150_p50_circle_SE_rounded = round2(std(SP_n150_p50_circle_vector)/sqrt(length(SP_n150_p50_circle_vector)),.01);

%have to exclude NaNs when reporting MCC.
MCC_n150_p50_circle_mean_rounded = round2(nanmean(MCC_n150_p50_circle_vector),.01);
MCC_n150_p50_circle_SE_rounded = round2(nanstd(MCC_n150_p50_circle_vector)/sqrt(length(MCC_n150_p50_circle_vector(~isnan(MCC_n150_p50_circle_vector)))),.01);


save('Bsplines_paper_Calcs_p50_n150_circle_SpikeSlab_1to100.mat',...
    'total_time_n150_p50_circle_mean_rounded',... 
'total_time_n150_p50_circle_SE_rounded',...
'L1_n150_p50_circle_mean_rounded',... 
'L1_n150_p50_circle_SE_rounded',...
   'SE_n150_p50_circle_mean_rounded',...
'SE_n150_p50_circle_SE_rounded', 'SP_n150_p50_circle_mean_rounded', 'SP_n150_p50_circle_SE_rounded', 'MCC_n150_p50_circle_mean_rounded',...
 'MCC_n150_p50_circle_SE_rounded');

