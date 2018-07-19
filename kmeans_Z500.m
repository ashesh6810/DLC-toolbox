% function [Sval] = kmeanscluster_PHZ(year1,year2,day1,day2,y1,y2,x1,x2,p1,p2,Swind,demean,manual)
clc;
clear all;
close all;
c=xmlread('input2.xml');
file_grid = ( c.getElementsByTagName( 'filename_grid' ).item( 0 ).getFirstChild.getNodeValue( ) );
filename_Z500= ( c.getElementsByTagName( 'filename_Z500' ).item( 0 ).getFirstChild.getNodeValue( ) );
load([char(file_grid) '/grid.mat'])
manual = 1;
ensembles=str2double( c.getElementsByTagName( 'ensembles' ).item( 0 ).getFirstChild.getNodeValue( ) );
cnt_ens=1; 
lat_north_index=str2double( c.getElementsByTagName( 'lat_north_index' ).item( 0 ).getFirstChild.getNodeValue( ) );
lat_south_index=str2double( c.getElementsByTagName( 'lat_south_index' ).item( 0 ).getFirstChild.getNodeValue( ) );
lon_west_index=str2double( c.getElementsByTagName( 'lon_west_index' ).item( 0 ).getFirstChild.getNodeValue( ) );
lon_east_index=str2double( c.getElementsByTagName( 'lon_east_index' ).item( 0 ).getFirstChild.getNodeValue( ) );
lat1=lat(lat_north_index+1:end);
[qx,qy]=meshgrid(lon(lon_west_index:lon_east_index),lat1(lat_south_index:end));



tic;
%% load Z500 and its anomalies
for m=1:ensembles
    load ([ char(filename_Z500) '/Z99daily_NA_M' num2str(cnt_ens) '.mat'])
    Zave=squeeze(mean(Z99NApattern(:,:,:,18:109),2));
    R{m}=Z99NApattern(:,:,:,18:109);
    for i=1:97
       anomalies(:,i,:,:)=squeeze(Z99NApattern(:,i,:,18:109))-Zave; 
    end
    M{m}=anomalies;
    cnt_ens=cnt_ens+1;

end

%% generate samples of vectorized Z500
count=1;
for m=1:ensembles
for i=61:86
    for k=1:92
       X(:,count) =reshape(M{m}(i,:,:,k),97*66,1);
       Psi(:,count)=reshape(R{m}(i,:,:,k),97*66,1);
count=count+1;
    end
end
end
%Do PCA
[EOFs,PCval]=EOFanalysis(X);


%if we want to manually select clusters from plot
if manual == 1
   % reply = input('How many clusters do you want?')   %enter number into command window and then hit 'Enter'
    %reply = input('How many clusters do you want?')   %enter number into command window and then hit 'Enter'
    %nC = reply;
   nC=str2double( c.getElementsByTagName( 'nC' ).item( 0 ).getFirstChild.getNodeValue( ) ); %%Change this to the number of classes you want   
    
    %reply = input('How many EOFs do you want?')   %enter number into command window and then hit 'Enter'
    %reply = input('How many EOFs do you want?')   %enter number into command window and then hit 'Enter'
    %nEOF = reply;
    nEOF = str2double( c.getElementsByTagName( 'nEOF' ).item( 0 ).getFirstChild.getNodeValue( ) ); %% change this to the number of EOFs you want to retain
    
%if we want to just select highest avg silhouette value between 4-20 clusters
elseif manual == 0
    EOFmax=30;
    Sval = zeros(4,EOFmax);
    for nEOF=5:EOFmax
        nEOF
        Xr = squeeze(EOFs(:,end-nEOF+1:end))'*X;
        Xtr=Xr';
        
        %calculate mean silhouette values for certain numbers of clusters
        silval = [];
        for r = 2:4
            [idx,C] = kmeans(Xtr,r,'replicates',100);
            S = silhouette(Xtr,idx);        %calculates silhouette values
            silval = [silval; r mean(S)];   %puts average silhouette value in table
            Sval(r,nEOF) = mean(S);
        end
    end
    h=figure(1)
    pcolor(1:4,1:EOFmax,Sval');colorbar
    im=frame2im(getframe(gca));
    imwrite(im,['silhouttevalues' num2str(ensembles) '.png'])
    disp('Chosen based on the silhouette values')
    [nC,nEOF] = find(max(max(Sval)) == Sval)
    close(h);
else
    disp('manual must be 0 or 1, no or yes')
end

Xr = squeeze(EOFs(:,end-nEOF+1:end))'*X;
Xtr=Xr';
sum(PCval(end-nEOF+1:end))*100.0/sum(PCval)

%kmeans replicated 1000 times (like Souri)
[idx, Cr] = kmeans(Xtr,nC,'replicates',1000);

Count(nC,1)=0;
for n=1:nC
    for d=1:length(idx)
        if(idx(d)==n)
            Count(n)=Count(n)+1;
        end
    end
end
[sum(Count) length(idx)]
h=figure(1)
silhouette(Xtr,idx);
saveas(h,['silhouttevalues' num2str(ensembles) '.png'])
close(h);
C = squeeze(EOFs(:,end-nEOF+1:end))*Cr';
C = C';

for j=1:nC
    count=1;
 for i=1:length(idx)

    if (idx(i)==j)
     cluster{j}(:,count)=Psi(:,i);
     count=count+1;
    end
 end
end
for j=1:4
 Z500_cn(:,j)=mean(cluster{j},2);
end
load coastlines

h=figure(4)
for j=1:nC
    ZZZ=reshape(Z500_cn(:,j),97,66);
    subplot(2,3,j)
    contourf(qx',qy',ZZZ,10,'LineColor','k','LineWidth',2);
    caxis([5400 5900])     
    hold on;

    plot(coastlon+360,coastlat,'Linewidth',1,'Color','k');
    xlim([195 315])
    ylim([25 100])
end
savefig(h,'Z500centers.fig')
close(h);

%plot code is generalized to stepping in data for lat/long
h=figure(2)
for n=1:nC
    subplot(ceil(nC/3),3,n)                 %makes subplots big enough
    Z=(reshape(C(n,1:size(C,2)),97,66));
    contourf(qx',qy',Z,10);hold on
    plot(coastlon+360,coastlat,'Linewidth',1,'Color','k');
    xlim([195 315])
    ylim([25 90])
    caxis([-120 120]);
    title(['Class' num2str(n) ' ' num2str(100*Count(n)/sum(Count))]);
    
end

%% save the cluster center images
savefig(h,'Z500anomalies.fig')
close(h)
hold off
%% plot the 1st 6 EOFs
h=figure(3)
for i=1:6
    k=i-1;
    subplot(3,2,i)
    contourf(qx',qy',reshape(EOFs(:,end-k),97,66),10,'LineColor','k','LineWidth',2);hold on;
    plot(coastlon+360,coastlat,'Linewidth',1,'Color','k');
    plot(coastlon+360,coastlat,'Linewidth',1,'Color','k');
    xlim([195 315])
    ylim([25 90])
    
end

savefig(h,['EOFS' num2str(nC)  '.fig'])
close(h)

toc

%% save the clustered data
save(['Clustereddatawithensemble' num2str(ensembles) 'withcluster' num2str(nC) '.mat'],'X','idx','Count','-v7.3')

% end
