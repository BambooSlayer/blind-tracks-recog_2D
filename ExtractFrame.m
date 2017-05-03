% function ExtractFrame(imagefile)
imagefile='MVI_8229.avi';
mov=aviread(imagefile);

%因为mov是个结构1*X的结构体，X就是帧的个数
%matlab是如何计算帧的个数，我不知道
%如果大家按照每秒25或30帧的个数截取，还是用软件的好
%size(mov,2)得出视频由多少帧组成
n=size(mov,2);

for i=1:n
    k=int2str(i);
    %截取某一帧
    F=mov(1,i);

    %可以查阅frame2im的使用方法
    [f,map]=frame2im(F);

    %存取路径
    k1=strcat(k,'.jpg');
    imwrite(f,k1);
end
