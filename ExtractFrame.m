% function ExtractFrame(imagefile)
imagefile='MVI_8229.avi';
mov=aviread(imagefile);

%��Ϊmov�Ǹ��ṹ1*X�Ľṹ�壬X����֡�ĸ���
%matlab����μ���֡�ĸ������Ҳ�֪��
%�����Ұ���ÿ��25��30֡�ĸ�����ȡ������������ĺ�
%size(mov,2)�ó���Ƶ�ɶ���֡���
n=size(mov,2);

for i=1:n
    k=int2str(i);
    %��ȡĳһ֡
    F=mov(1,i);

    %���Բ���frame2im��ʹ�÷���
    [f,map]=frame2im(F);

    %��ȡ·��
    k1=strcat(k,'.jpg');
    imwrite(f,k1);
end
