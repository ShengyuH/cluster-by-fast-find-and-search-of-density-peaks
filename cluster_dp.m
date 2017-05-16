function [cluster_lables, center_idxs,rho,delta] = cluster_dp(dist, para,data,num_MSE)
%% Input and Output
% INPUT :
% dist : A nCases*nCases matrix. each dist(i, j) represents the distance
%        between the i-th datapoint and j-th datapoint.
% para : options
%        percent : The parameter to determine dc. 1.0 to 2.0 can often yield good performance.
%        method  : alternative ways to compute density. 'gaussian' or
%                  'cut_off'. 'gaussian' often yields better performance.
%num_MSE:used in AutoSet cluster centroids
%data:   original data used in building distance matrix
% OUTPUT :
% cluster_labels : a nCases vector contains the cluster labls. Lable equals to 0 means it's in the halo region
% center_idxs    : a nCluster vector contains the center indexes.


%% Estimate dc
disp('Estimating dc...');
percent = para.percent;
N = size(dist,1);
position = round(N*(N-1)*percent/100);
tri_u = triu(dist,1);
sda = sort(tri_u(tri_u~=0));
dc = sda(position);
fprintf('Computing dc with original method: %12.6f\n', dc);
clear sda; clear tri_u;

%% new method to estimate dc
% disp('Estimating dc...');
% STD=std2(data);
% dc=(4/(3*size(dist,1)))^(1/5)*STD;
% fprintf('Computing dc with new method: %12.6f\n', dc);

%% Compute rho(density)
disp('Computing rho...');
switch para.method
    % Gaussian kernel
    case 'gaussian'
        rho = sum(exp(-(dist./dc).^2),2)-1;
        % "Cut off" kernel
    case 'cut_off'
        rho = sum((dist-dc)<0, 2);
end
[~,ordrho]=sort(rho,'descend');

rho(rho(:)==0)=0.0000001;
%% Compute delta
disp('Computing delta...');

delta = zeros(size(rho));
nneigh = zeros(size(rho));

delta(ordrho(1)) = -1;
nneigh(ordrho(1)) = 0;
for i = 2:size(dist,1)
    range = ordrho(1:i-1);
    [delta(ordrho(i)), tmp_idx] = min(dist(ordrho(i),range));
    nneigh(ordrho(i)) = range(tmp_idx); 
end
delta(ordrho(1)) = max(delta(:));
delta=delta+1;

%% normalize delta and rho
%rho=mapminmax(rho',0.0000001,1)';
%delta=mapminmax(delta',0.000001,1)';

%% %fit (rho-delta)  y=a*x^b;
[paras,S]=polyfit(log(delta),log(rho),1);
estrho=exp(polyval(paras,log(delta)));
MSE_rho=sqrt((S.normr)^2/S.df);
disp(['The Mean Square Residuals is:',num2str(MSE_rho)]);

[paras,S]=polyfit(log(rho),log(delta),1);
estDelta=exp(polyval(paras,log(rho)));
MSE_delta=sqrt((S.normr)^2/S.df);
disp(['The Mean Square Residuals is:',num2str(MSE_delta)]);


%% Decision graph, choose min_rho and min_delta
figure(1)
plot(rho(:),delta(:),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
C=sortrows([rho,estDelta],1);
plot(C(:,1),C(:,2),'r',...
    'LineWidth',2);
hold on
C=sortrows([rho,estDelta+MSE_delta*num_MSE],1);
plot(C(:,1),C(:,2),'g',...
    'LineWidth',2);
hold on
C=sortrows([estrho,delta],1);
plot(C(:,1),C(:,2),'r--',...
    'LineWidth',2);
hold on
C=sortrows([estrho+MSE_rho*num_MSE,delta],1);
plot(C(:,1),C(:,2),'g--',...
    'LineWidth',2);
    
%% Rectset rho_min,delta_min
% figure(1)
% rect = getrect(1); 
% rho_min=rect(1);
% delta_min=rect(2);

%% Autoset rho_min,delta_min
delta_min=estDelta+MSE_delta*num_MSE;
rho_min=estrho+MSE_rho*num_MSE;
% c=[min(delta):0.01:max(delta)];
% r=rho_min*ones(1,size(c,2));
% hold on
% plot(r,c,'c',...
%     'linewidth',2);


%% Find cluster centers
center_idxs=intersect(find(rho>rho_min),find(delta>delta_min));
while(true)   
%% Decision graph, choose min_rho and min_delta
disp('Finding cluster centers...');
disp([num2str(length(center_idxs)),' cluster centers found...']);
close all
figure(1)
plot(rho(:),delta(:),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
C=sortrows([rho,estDelta],1);
plot(C(:,1),C(:,2),'r',...
    'LineWidth',2);
hold on
C=sortrows([rho,estDelta+MSE_delta*num_MSE],1);
plot(C(:,1),C(:,2),'g',...
    'LineWidth',2);
hold on
C=sortrows([estrho,delta],1);
plot(C(:,1),C(:,2),'r--',...
    'LineWidth',2);
hold on
C=sortrows([estrho+MSE_rho*num_MSE,delta],1);
plot(C(:,1),C(:,2),'g--',...
    'LineWidth',2);
    
%% show cluster centers
NCLUST=length(center_idxs);
cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   hold on
   plot(rho(center_idxs(i)),delta(center_idxs(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
%% Assignment
% raw assignment
    disp('Assigning data-points into clusters...');
    cluster_lables = -1*ones(size(dist,1),1);
    for i = 1:length(center_idxs)
        cluster_lables(center_idxs(i)) = i;
    end
    for i=1:length(cluster_lables)
        if (cluster_lables(ordrho(i))==-1)
            cluster_lables(ordrho(i)) = cluster_lables(nneigh(ordrho(i)));
        end
    end

%% rule out wrong cluster centroids
    fretable=tabulate(cluster_lables(:,1)); 
    wrongcenters=find(fretable(:,3)<2);
    if(size(wrongcenters,1)==0)
        box off
        break
    else
        orgindex=1:NCLUST;
        newindex=setdiff(orgindex,wrongcenters);
        center_idxs=center_idxs(newindex);
    end
    box off
    c=1;
end

raw_cluster_lables = cluster_lables;
%% find and cut off halo region
% disp('Cut off halo regions...');
% for i = 1:length(center_idxs)
%     tmp_idx = find(raw_cluster_lables==i);
%     tmp_dist = dist(tmp_idx,:);
%     tmp_dist(:,tmp_idx) = max(dist(:));
%     tmp_rho = rho(tmp_idx);
%     tmp_lables = raw_cluster_lables(tmp_idx);
%     tmp_border = find(sum(tmp_dist<dc,2)>0);
%     if ~isempty(tmp_border)
%         rho_b = max(tmp_rho(tmp_border));
%         halo_idx = rho(tmp_idx) < rho_b;
%         tmp_lables(halo_idx) = 0;
%         %lable equals to 0 means it's in the halo region
%         cluster_lables(tmp_idx) = tmp_lables;
%     end
% end

end