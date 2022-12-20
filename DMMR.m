function [cover_lsb,change,rho] = DMMR(cover_Path,rho1_P,rho1_M,tab_m)
%% �������ƻ��Ƽ�������Ԫ�ص��޸�ʧ����޸ľ���

cover_spa = imread(cover_Path);  
bits = 8;            %��������ֵת��Ϊ������DCTϵ��Ҫ������ֵ����ƽ�ƣ�ʹ����ֵ�ķֲ���0�Գ�                       
cover_spa = double(cover_spa) - 2^(round(bits) - 1);
[xm,xn] = size(cover_spa);
fun = @(x)dct2(x.data);
cover_DCT = blockproc(cover_spa,[8 8],fun);
m_block = floor(xm/8);
n_block = floor(xn/8);

%% �������ƻ��Ƽ������
G = 1;
usable_DCT_num = 9;
n_lsb = 0;
cover_round = zeros(1,m_block*n_block* usable_DCT_num); 
cover_lsb = zeros(1,m_block*n_block* usable_DCT_num);
change = zeros(1,m_block*n_block* usable_DCT_num); %Ƕ��ʱ����DCTϵ�� +- �ı���
rho = zeros(1,m_block*n_block* usable_DCT_num);   %Ƕ��ʱDCTϵ�� +-1 ��ʧ��
rho1_P_deq  = zeros(1,m_block*n_block* usable_DCT_num); %������DCTϵ�� +-1 ��ʧ��
rho1_M_deq = zeros(1,m_block*n_block* usable_DCT_num);

for bm = 1:m_block
    for bn = 1:n_block
        for i = 1:8
            for j = 1:8
                if (i+j==5)||(i+j==6)%��Ƶ9��DCTϵ��
                    n_lsb = n_lsb + 1;
                    rho1_P_deq(n_lsb) = rho1_P( (bm-1)*8+i,(bn-1)*8+j ) / double( tab_m(i,j) );  %������DCTϵ��+-1ʧ��
                    rho1_M_deq(n_lsb) = rho1_M((bm-1)*8+i,(bn-1)*8+j)/double(tab_m(i,j));
                    yd = cover_DCT((bm-1)*8+i,(bn-1)*8+j);  %��������DCTϵ��
                    tab_q = double(tab_m(i,j))/G;  %��������
                    cover_round(n_lsb)=round(yd/tab_q); 
                    dnum1 = round(yd/tab_q);
                     if mod(dnum1,2)==0
                        cover_lsb(n_lsb)=0;
                        dnum2 = floor(yd/tab_q);
                        if mod(dnum2,2)==1
                            change(n_lsb) = -1*(yd-dnum2*tab_q);%��ȥ��ô��
                            rho(n_lsb) = -1*change(n_lsb)*rho1_M_deq(n_lsb);
                        else
                            change(n_lsb) = ((dnum2+1)*tab_q-yd);%������ô��
                            rho(n_lsb) = change(n_lsb)*rho1_P_deq(n_lsb);
                        end                        
                    else
                        cover_lsb(n_lsb)=1;
                        dnum2 = floor(yd/tab_q);
                        if mod(dnum2,2)==1
                            change(n_lsb) = (dnum2+1)*tab_q-yd;%������ô��
                            rho(n_lsb) = change(n_lsb)*rho1_P_deq(n_lsb);
                        else
                            change(n_lsb) = -1*(yd-dnum2*tab_q);%������ô��
                            rho(n_lsb) = -1*change(n_lsb)*rho1_M_deq(n_lsb);
                        end
                    end  
                    
                end
            end
        end
    end
end

end