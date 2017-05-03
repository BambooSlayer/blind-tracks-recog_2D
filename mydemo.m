%% ��ƵԤ���� 
hsrc = vision.VideoFileReader('MVI_8229.avi', ...
                                  'VideoOutputDataType', 'uint8');
    hfg = vision.ForegroundDetector(...
        'NumTrainingFrames', 150, ... % 5 because of short video
        'InitialVariance', (30)^2,...
        'NumGaussians',4,...
        'MinimumBackgroundRatio',0.7); % initial standard deviation of 30
    hfg1 = vision.ForegroundDetector(...
        'NumTrainingFrames', 150, ... % 5 because of short video
        'InitialVariance', (30)^2,...
        'NumGaussians',4,...
        'MinimumBackgroundRatio',0.7); % initial standard deviation of 30
    hblob = vision.BlobAnalysis(...
        'CentroidOutputPort', false, 'AreaOutputPort', false, ...%������������꣬��������
        'BoundingBoxOutputPort', true, 'MinimumBlobArea', 500);
    hsi = vision.ShapeInserter('BorderColor','White');
  ti=20;
  n=1;
%     hsnk = vision.VideoPlayer();
    while ti~=1%~isDone(hsrc)
      step(hsrc);%�����Ƶ�е���һ֡
       step(hsrc);%�����Ƶ�е���һ֡
      frame  = step(hsrc);%�����Ƶ�е���һ֡
     %%
    ti=ti-1;
    if(ti==0)%����ģ��       
        reset(hfg1);
        ti=20;
    end

     %% 
      fgMask = step(hfg, frame); %���ǰ����mask����ʽ�ǲ�����
      fgMask1 = step(hfg1, frame);
     % bbox= step(hblob, fgMask);
      mask=fgMask.*ones(480,640)*255;
      mask1=fgMask1.*ones(480,640)*255;
      maskr=frame;
      maskr(:,:,1)=mask;
      maskr(:,:,2)=mask1;
      maskr(:,:,3)=zeros(480,640);
       if(ti==1)%
        PR(n,:,:)=mask;
        PG(n,:,:)=mask1;
%         n=n+1;
        break
%      out   = frame+maskr; % draw
%       step(hsnk, out); % view results in the video player
    end
    
    end
%     release(hsnk);
    release(hsrc);
%% ��Ƶ������� 
PG2=zeros(480,640);
PG2(:,:)=PG(1,:,:).*PR(n,:,:);%��϶���ѧϰ�볤��ѧϰ��õ�ǰ��
x=-90:640;
%% ¼��ͼ����ʾ
f=PG2;
figure(1);
subplot(2,2,1);imshow(frame+maskr);title('ԭͼ');

%% ��ȡͼ���Ե
[m,n]=size(f);%�õ�ͼ���������m������n
for i=3:m-2
    for j=3:n-2%��������ϴ����Դ�ͼ��3,3����ʼ���ڣ�m-2,n-2������ 
        l(i,j)=-f(i-2,j)-f(i-1,j-1)-2*f(i-1,j)-f(i-1,j+1)-f(i,j-2)-2*f(i,j-1)+16*f(i,j)-2*f(i,j+1)-f(i,j+2)-f(i+1,j-1)-2*f(i+1,j)-f(i+1,j+1)-f(i+2,j);%LoG����
    end
end
subplot(2,2,2);imshow(l);title('LoG������ȡͼ���Ե');

%%  �˲�
[m,n]=size(l);
for i=2:m-1    
 for j=2:n-1   
    y(i,j)=l(i-1,j-1)+l(i-1,j)+l(i-1,j+1)+l(i,j-1)+l(i,j)+l(i,j+1)+l(i+1,j-1)+l(i+1,j)+l(i+1,j+1);
    y(i,j)=y(i,j)/9;  %LoG������ȡ��Ե�󣬶Խ�����о�ֵ�˲���ȥ��������Ϊ��һ��hough�任��ȡֱ����׼��
 end
end
subplot(2,2,3);imshow(y);title('��ֵ�˲��������')

%% ��ֵ��
q=im2uint8(y);
[m,n]=size(q);
for i=1:m    
 for j=1:n  
    if q(i,j)>80;   %���ö�ֵ������ֵΪ80
        q(i,j)=255; %��ͼ����ж�ֵ������ʹͼ���Ե����ͻ������
    else
        q(i,j)=0; 
    end
 end
end
subplot(2,2,4);imshow(q);title('��ֵ�������');

%% ���ֱ��
 
%Hough�任���ֱ�ߣ�ʹ�ã�a��p�������ռ䣬a��[0,180],p��[0,2d]
a=180; %�Ƕȵ�ֵΪ0��180��
d=round(sqrt(m^2+n^2)); %ͼ��Խ��߳���Ϊp�����ֵ
s=zeros(a,2*d); %�洢ÿ��(a,p)����
[aa,bb]=find(q==255);
z=zeros(1,2);   %bu��Ԫ���洢ÿ�������ĵ������
zz=zeros(1,2);
zzz=zeros(1,2);%�������ļ�
ll=length(aa);%ll��
for loop=1:ll
        i=aa(loop);
     j=bb(loop);
            for k=1:a
                p = round(i*cos(pi*k/180)+j*sin(pi*k/180));%��ÿ�����1��180�ȱ���һ�飬ȡ�þ����õ������ֱ�ߵ�pֵ(����)��ȡ����
                if(p > 0)%��p����0���򽫵�洢�ڣ�d��2d���ռ�
                    s(k,d+p)=s(k,d+p)+1;%��a��p����Ӧ���ۼ�����Ԫ��һ
                    if(s(k,d+p) == 125)
                       z=[z;k,p];
                    end
                else
                    ap=abs(p)+1;%��pС��0���򽫵�洢�ڣ�0��d���ռ�
                    s(k,ap)=s(k,ap)+1;
                   if(s(k,ap) ==125)
                    z=[z;k-180,ap]; 
                    end 
                end
            end
       
end


%% ��ʾЧ��
figure(22);                                                        
plot(0,0);
axis([0 640 0 480]);
title('ֱ�߼��');
hold;

% xx=z(:,2)./sin((pi.*z(:,1)./180));
yy=z(:,2)./sin(pi.*z(:,1)./180);
k=-cot(pi.*z(:,1)./180);

a=find(z(:,1)~=0);
for i=2:(1+length(a))    
 y=k(i).*x+yy(i);
 plot(y,480-x,'-r','LineWidth',1); %%%%%%%%%%%%%%%55
  D(i)=y(1);
end

figure(201)
for i=-800:800
    DD(i+801)=sum(D(D>((i-1)*10+1)&D<i*10));
end

bar(-800:800,DD)
title('����ƽ���ϵĽؾ�ͳ��');
xlabel('��ƽ�߽ؾ�')
ylabel('ֱ����Ŀͳ��')

%%
figure(202)
plot(z(:,1),z(:,2),'r+')
hold on
i=find(DD==max(DD(801:801+70)))-801;%
zz=[z((D>((i-1)*10+1)&D<i*10),1),z((D>((i-1)*10+1)&D<i*10),2)];
plot(zz(:,1),zz(:,2),'go')%Ȧ����ʧ����ǰ���ĸ��������ĵ�
axis([-30 180 0 480]);
title('��������ֱ�ߴ�����ֵ');
xlabel('�Ƕ�')
xlabel('����')
%% ʶ����
figure(203)
a=find(zz(:,1)~=0);
xx=zz(:,2)./sin((pi.*zz(:,1)./180));
yy=zz(:,2)./sin(pi.*zz(:,1)./180);
k=-cot(pi.*zz(:,1)./180);
title('ɸѡ���ֱ��');
imshow(frame);hold on;
for i=1:(length(a))    
 y=k(i).*x+yy(i);
 plot(y,x,'g')
end
axis([0 640 0 480]);
%%