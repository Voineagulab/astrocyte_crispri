clear;
%%%%%%%%%%%%%%%%%%%%%% CONFIG %%%%%%%%%%%%%%%%%%%%%%
TRD_number_of_read=2;
TRD_pvalue=3; %3= ln plvaue (0.05); 4.6=ln p-value (0.01); 6.9 = ln p-value (0.001)
binsize='../bin_sizes/5000/';
filename='cis_interactions.txt';
%%%%%%%%%%%%%%%%%%%%%% CONFIG %%%%%%%%%%%%%%%%%%%%%%

%%% import interactions from MAXHiC %%%
fprintf('importing interactions ...\n');
formatstring='%d %d %d %f %f %*[^\n]';
fid=fopen(filename, 'r');
sample_interaction=textscan(fid, formatstring, 'delimiter', ',','headerLines', 1);
fclose(fid);

%%% filter intractions based on nubmer of reads %%%
%indx1=sample_interaction{1,3}>= TRD_number_of_read;
%indx2=sample_interaction{1,4}>= TRD_pvalue;
%indx3=abs(sample_interaction{1,2}-sample_interaction{1,1})>0;
%indx= indx1 & indx2 & indx3;
%sig_interaction = cellfun(@(x) x(indx), sample_interaction, 'UniformOutput', false);
sig_interaction=sample_interaction;

%%% import Astrocytes enhancer midpoint gene %%%
fprintf('importing enhancer data...\n');
inputfile='Astrocytes_HG38_Enh_Midpoint_Gene_TSS_FID.txt';
fid=fopen(inputfile, 'r');
formatstring='%s %s %d %d %s %*d %*d %*d %*d %d %d %*[^\n]';
%Final_matrix=textscan(fid, formatstring, 1, 'delimiter', '\t');
Enhancer_info=textscan(fid, formatstring, 'delimiter', '\t', 'HeaderLines', 1);
fclose(fid);
Final_matrix=cell(1,6);
Final_matrix(1,1)=cellstr('Enh.chr'); Final_matrix(1,2)=cellstr('Enh.Midpoint'); Final_matrix(1,3)=cellstr('Gene.TSS');
Final_matrix(1,4)=cellstr('Pair'); Final_matrix(1,5)=cellstr('#HiC.interaction'); Final_matrix(1,6)=cellstr('P-value');

%%% Count number of interactions between enhancer and the gene TSS %%%
number_of_enhancers=length(Enhancer_info{1,1});
tic;
for i= 1 : number_of_enhancers
    fprintf('region %1.0f of %1.0f ...\n', i, number_of_enhancers);
    enhancer_point_FID=Enhancer_info{1,6}(i);
    gene_tss_FID=Enhancer_info{1,7}(i);
    fragmentID_left=sig_interaction{1,1};
    fragmentID_right=sig_interaction{1,2};
    indx1=find(sig_interaction{1,1}==enhancer_point_FID & sig_interaction{1,2}==gene_tss_FID);
    indx2=find(sig_interaction{1,1}==gene_tss_FID & sig_interaction{1,2}==enhancer_point_FID);
    
    Final_matrix(i+1,1)=Enhancer_info{1,2}(i); Final_matrix(i+1,2)=num2cell(Enhancer_info{1,3}(i));
    Final_matrix(i+1,3)=num2cell(Enhancer_info{1,4}(i)); Final_matrix(i+1,4)=num2cell(Enhancer_info{1,5}(i));
    Final_matrix(i+1,5)=num2cell(sum(sig_interaction{1,3}(indx1))+ sum(sig_interaction{1,3}(indx2)));
    if indx1
        Final_matrix(i+1,6)=num2cell(max(sig_interaction{1,4}(indx1)));
    elseif indx2
        Final_matrix(i+1,6)=num2cell(max(sig_interaction{1,4}(indx2)));
    else
        Final_matrix(i+1,6)=num2cell(0);
    end 
end
toc
