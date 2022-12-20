function [sh_I] = Image_Shuffle(I,SH)
% ����˵������ͼ��I���ݻ�ϴ����SH���л�ϴ
% ���룺I��ԭʼͼ�����,SH����ϴ���У�
% �����sh_I����ϴ���ͼ�����

[m,n] = size(I);
Shuffle = reshape(SH,m,n); %���ɾ���,��������
x_I = zeros(1,numel(Shuffle)); %����һά��ϴ����
for j=1:n
    for i=1:m
        x = Shuffle(i,j); %��i,j�����ҵ�λ������Ϊx
        x_I(x) = I(i,j);
    end
end
sh_I = reshape(x_I,m,n);
end