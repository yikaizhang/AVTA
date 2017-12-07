%********This program runs the swimmer image NMF experiment*********
%Y.mat is the original data set which should be included in current path.
% Compare AVTA and Fast anchor word on NMF problem. The goal is to recover
% under lying poses of different limbs. Data are corrupted by spurious
% action.
%*******************************************************************


%% load data
swimraw=load('Y.mat');
Y_data=swimraw.Y;


%% Generate corrupted data
[pi,pj,np]=size(Y_data);
YY=Y_data;
noise_v=3;
action=YY(:,:,:);
scale=0.5;%%Scale of noise
for i=1:np
    YY(:,:,i)=Y_data(:,:,i)+randblock(action(:,:,randi(256)),[32,8])*scale+rand(32);
    %%% The data is corrupted by a random block of a image and some iid
    %%% uniform random noise.
end

%% Translation: Image -> Word-Document vector
word_doc_mat=zeros(pi*pj,np);
for i=1:np
    Yi=YY(:,:,i);
    word_doc_mat(:,i)=Yi(:);
end

%% Compute underlying topic
K=16;

[aa_A]=recover_l2(word_doc_mat,K);
[at_A]=recover_l2_avta(word_doc_mat,K);





%% Translation: Doc-Word vector -> Image
image_at=zeros(32,32,16);
image_aa=zeros(32,32,16);
image_tsvd=zeros(32,32,16);
image_tvt=zeros(32,32,16);
for i=1:16
    image_at_i=reshape(at_A(:,i),32,32)*2430;
    image_aa_i=reshape(at_A(:,i),32,32)*2430;
    image_at(:,:,i)=image_at_i;
    image_aa(:,:,i)=image_aa_i;
end
% 
% 
% 
% %% Output result One 2x16 figure
% hold on;
% for j=1:16
%     subplot(2,16,j);
%     im=image_aa(:,:,j);
%     im(image_aa(:,:,j)<40)=0;%%Threshold output
%     image(im);
% 
%     set(gca,'xtick',[],'ytick',[])
% end
% 
% 
% 
% for j=1:16
%     subplot(2,16,16+j);
%     im=image_at(:,:,j);
%     im(image_at(:,:,j)<40)=0;
%     image(im);
%     set(gca,'xtick',[],'ytick',[])
% end
% 
% hold off;
% 
% 

% %% Output result Two 4x4 figures
% 
% hold on;
% suptitle('Separable NMF +RecoverL2')
% for j=1:16
%     subplot(4,4,j);
%     
%     im=image_aa(:,:,j);
%     im(image_aa(:,:,j)<40)=0;
%     image(im);
% 
%     set(gca,'xtick',[],'ytick',[])
% end
% 
% legend('show','Location','east')
% hold off;
% hold on;
% suptitle('AVTA+RecoverL2')
% for j=1:16
%     subplot(4,4,j);
%     im=image_at(:,:,j);
%     im(image_at(:,:,j)<40)=0;
%     image(im);
%     set(gca,'xtick',[],'ytick',[])
% end
% 
% hold off;
% 
% legend('show','Location','east')


%% Plot swimmer dataset and topics: Press anykey to continue
hold on;
suptitle('Separable NMF +RecoverL2')
for j=1:16
    subplot(4,4,j);
    
    im=image_aa(:,:,j);
    im(image_aa(:,:,j)<40)=0;
    image(im);

    set(gca,'xtick',[],'ytick',[])
end

legend('show','Location','east')
hold off;
pause();
hold on;
suptitle('AVTA+RecoverL2')
for j=1:16
    subplot(4,4,j);
    im=image_at(:,:,j);
    im(image_at(:,:,j)<40)=0;
    image(im);
    set(gca,'xtick',[],'ytick',[])
end

hold off;

legend('show','Location','east')


pause();

suptitle('swimmer image');
pause();
hold on;
for i=5:8
    for j=1:16
        subplot(8,8,16*(i-5)+j);
        image(YY(:,:,16*(i-1)+j));
        set(gca,'xtick',[],'ytick',[])
    end
end
hold off;
suptitle('swimmer image');
pause();
hold on;
for i=9:12
    for j=1:16
        subplot(8,8,16*(i-9)+j);
        image(YY(:,:,16*(i-1)+j));
        set(gca,'xtick',[],'ytick',[])
    end
end
hold off;
suptitle('swimmer image');
pause();
hold on;
for i=13:16
    for j=1:16
        subplot(8,8,16*(i-13)+j);
        image(YY(:,:,16*(i-1)+j));
        set(gca,'xtick',[],'ytick',[])
    end
end
hold off;
suptitle('swimmer image');
pause();


