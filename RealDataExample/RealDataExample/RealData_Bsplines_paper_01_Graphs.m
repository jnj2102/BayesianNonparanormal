%Code to create the estimated graphs for the Real Data application
%
%Author: Jami Jackson Mulgrave

load('RealData_Bsplines_paper_nonparanormal_stars_ric_ebic_01.mat');
load('Bsplines_paper_realdata_01.mat');

%Bayesian nonparanormal graph
edge_matrix_realdata_finalanalysis(logical(eye(size(edge_matrix_realdata_finalanalysis)))) = 0;

RealData_bayes = graph(edge_matrix_realdata_finalanalysis, 'upper');



H = plot(RealData_bayes); %I can remove the axes with the plotting tool.
%set(gca,'xtick',[],'ytick',[])

H.NodeLabel = {};  %this removes the node labels altogether.
set(gca, 'XTick',[])
set(gca, 'YTick',[])  %this removes the tick marks - now it rivals graphviz4matlab.
set(gca,'Visible','off') %remove axes border

saveas(H, 'RealData_bayes_01_nonames.jpg')

%Read in the frequentist STARS nonpararnormal model

edge_matrix_nonparanormal_stars(logical(eye(size(edge_matrix_nonparanormal_stars)))) = 0;

RealData_nonparanormal_stars = graph(edge_matrix_nonparanormal_stars, 'upper');

F = plot(RealData_nonparanormal_stars);

F.NodeLabel = {};  %this removes the node labels altogether.
set(gca, 'XTick',[])
set(gca, 'YTick',[])  %this removes the tick marks - now it rivals graphviz4matlab.
set(gca,'Visible','off')

saveas(F, 'RealData_nonparanormal_stars_01_nonames.jpg')

%Read in the frequentist STARS alternative (threshold = .05) nonpararnormal model

edge_matrix_nonparanormal_stars_alternative(logical(eye(size(edge_matrix_nonparanormal_stars_alternative)))) = 0;

RealData_nonparanormal_stars_alternative = graph(edge_matrix_nonparanormal_stars_alternative, node_names, 'upper');

K = plot(RealData_nonparanormal_stars_alternative);

K.NodeLabel = {};  %this removes the node labels altogether.
set(gca, 'XTick',[])
set(gca, 'YTick',[])  %this removes the tick marks - now it rivals graphviz4matlab.
set(gca,'Visible','off')


saveas(K, 'RealData_nonparanormal_stars_alternative_01.jpg')

%Read in the frequentist RIC nonpararnormal model

edge_matrix_nonparanormal_ric(logical(eye(size(edge_matrix_nonparanormal_ric)))) = 0;

RealData_nonparanormal_ric = graph(edge_matrix_nonparanormal_ric,  'upper');

G = plot(RealData_nonparanormal_ric);

G.NodeLabel = {};  %this removes the node labels altogether.
set(gca, 'XTick',[])
set(gca, 'YTick',[])  %this removes the tick marks - now it rivals graphviz4matlab.
set(gca,'Visible','off')

saveas(G, 'RealData_nonparanormal_ric_01_nonames.jpg')

%Read in the frequentist EBIC nonpararnormal model

edge_matrix_nonparanormal_ebic(logical(eye(size(edge_matrix_nonparanormal_ebic)))) = 0;

RealData_nonparanormal_ebic = graph(edge_matrix_nonparanormal_ebic, 'upper');

I = plot(RealData_nonparanormal_ebic);

I.NodeLabel = {};  %this removes the node labels altogether.
set(gca, 'XTick',[])
set(gca, 'YTick',[])  %this removes the tick marks - now it rivals graphviz4matlab.
set(gca,'Visible','off')


saveas(I, 'RealData_nonparanormal_ebic_01_nonames.jpg')

%How many edges do we get for each?

[p, ~] = size(edge_matrix_realdata_finalanalysis);

indmx = reshape(1:p^2,p,p); 
  upperind = indmx(triu(indmx,1)>0);  %do not include the diagonal
 
  sum_edges_bayes =  sum(edge_matrix_realdata_finalanalysis(upperind) == 1);

    sum_edges_stars =  sum(edge_matrix_nonparanormal_stars(upperind) == 1);

      sum_edges_stars_alternative =  sum(edge_matrix_nonparanormal_stars_alternative(upperind) == 1);

        sum_edges_ric =  sum(edge_matrix_nonparanormal_ric(upperind) == 1);

          sum_edges_ebic =  sum(edge_matrix_nonparanormal_ebic(upperind) == 1);
