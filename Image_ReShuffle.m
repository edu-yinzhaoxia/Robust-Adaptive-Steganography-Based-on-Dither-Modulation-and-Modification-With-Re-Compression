function [re_I] = Image_ReShuffle(sh_I,SH)
% ����˵������ͼ��I���ݻ�ϴ����SH���л�ϴ
% ���룺sh_I����ϴ��ͼ�����,SH����ϴ���У�
% �����re_I���ָ����ͼ�����

[m,n] = size(sh_I);
num = numel(SH);
I1 = reshape(sh_I,1,m*n); %��sh_Iת����һά����
x_I = zeros(1,m*n); %����һά�ָ�����
for i=1:num
    x = SH(i); %���ҵ�����
    x_I(i) = I1(x);
end
re_I = reshape(x_I,m,n);
end