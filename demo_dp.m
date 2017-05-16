%% demo of clustering by finding density peaks
clear all
close all
clc

%% simulate data
% mu1=[10 10 10];  
% S1=[2 0 0;0 2 0;0 0 2];  
% data1=mvnrnd(mu1,S1,60);   
% 
% mu2=[30 30  10];
% S2=[2 0 0;0 2 0;0 0 2];
% data2=mvnrnd(mu2,S2,100);
% 
% mu3=[10 30 30];
% S3=[2 0 0;0 2.1 0;0 0 2.2];
% data3=mvnrnd(mu3,S3,200);
% 
% mu4=[30 10 30];
% S4=[2 0 0;0 2.2 0;0 0 2.1];
% data4=mvnrnd(mu4,S4,200);
% 
% data5=2*rand(3,3);
% 
% data=[data1;data2;data3;data4;data5];
% [h,w]=size(data);
% for i=1:w
%     data(:,i)=mapminmax(data(:,i)',0,1)';
% end

%% data from paper
% data=load('fig2_panelB.dat');
% [h,w]=size(data);
% for i=1:w
%     data(:,i)=mapminmax(data(:,i)',0,1)';
% end

%% read UCI dataset
xx=load('S3.dat');
data=xx(:,2:3);
[h,w]=size(data);
for i=1:w
   data(:,i)=mapminmax(data(:,i)',0,1)'; 
end

%% read normalImage or ENVIImage
%filename='file.tif';

%normal Image
%I=single(imread(filename));

%ENVIImageFile
% I=read_ENVIimagefile;
% 
% %histogram equalization
% [h,w,dim]=size(I);
% M=h*w;
% data=reshape(I,h*w,dim);
% % for i=1:dim
%    data(:,i)=mapminmax(data(:,i)',0,1)'; 
%    img=single(HistogramEqualization(reshape(data(:,i),h,w)));
%    data(:,i)=mapminmax(reshape(img,h*w,1)',0,1)';
% end

%% dimension reduction
% [coeff,score,latent]=pca(data);
% latent=latent./sum(latent);
% sum=0;
% i=0;
% while (sum<0.95)
%     i=i+1;
%     sum=sum+latent(i);
% end
% data=score(:,1:i);
% dim=size(data,2);
% for i=1:dim
%    data(:,i)=mapminmax(data(:,i)',0,1)'; 
% end

%% plot dist
% figure(2)
% D=pdist(data);
% c=downsample(D,12000);
% k=1:size(c,2);
% plot(k,c,'*');

%% Clustering
dist=pdist2(data,data);
para.method = 'gaussian';
para.percent = 2;
num_MSE=2.5;
[cluster_lables, center_idxs,rho,delta] = cluster_dp(dist, para,data,num_MSE);
clear dist

%cluster_lables=kmeans(data,5);
%show the result
% figure(2)
% % %brighten(-0.85)
% cmap = colormap;
% %nclust=5;
% nclust = length(center_idxs);
% for i = 1:nclust
%     tmp_data = data(cluster_lables==i,:);
%     ic = int8((i*64.)/(nclust*1.));
%     col = cmap(ic,:);
% %     plot3(tmp_data(:,1),tmp_data(:,2),tmp_data(:,3),...
% %          'o','MarkerEdgeColor','k','MarkerFaceColor',col); 
%     plot(tmp_data(:,1),tmp_data(:,2),...
%         'o','MarkerSize',2,'MarkerFaceColor',col,'MarkerEdgeColor',col);
%     hold on;
% end
% tmp_data = data(cluster_lables==0,:);
% %  plot3(tmp_data(:,1),tmp_data(:,2),tmp_data(:,3),...
% %          'o','MarkerEdgeColor','k','MarkerFaceColor','k'); 
%      plot(tmp_data(:,1),tmp_data(:,2),...
%         'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
%     t=0:0:0;
%     set(gca,'xtick',t)
%     set(gca,'ytick',t)
% static  
% cluster_1=find(cluster_lables==1);
% cluster_2=find(cluster_lables==2);
% cluster_3=find(cluster_lables==3);

%show image
% figure(3)
% I1=reshape(cluster_lables,h,w);
% imagesc(I1);
% colorbar
% t=0:0:0;
% set(gca,'xtick',t)
% set(gca,'ytick',t)
t=tabulate(cluster_lables);

    
    
