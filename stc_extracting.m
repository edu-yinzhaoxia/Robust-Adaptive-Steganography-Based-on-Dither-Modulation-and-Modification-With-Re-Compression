function [stc_decoded_msg] = stc_extracting(afterchannel_stego_Path,stc_n_msg_bits,stc_n_msg_bits_check,tab_m,n,k,m,l_check_check)
%% ��ȡ��������Ԫ��
    bits = 8;
    cover_spa = imread(afterchannel_stego_Path);
    cover_spa = double(cover_spa) - 2^(round(bits)-1);
    [xm,xn] = size(cover_spa);
    t = dctmtx(8);
    fun = @(xl) (t*xl*(t'));
    cover_DCT = blkproc(cover_spa,[8 8],fun);%�ֿ�DCT�任
    m_block = floor(xm/8);
    n_block = floor(xn/8);
    G = 1;
    n_msg = 0;
    code_n = m_block*n_block*9;  %��Ƶ9��DCTϵ��
    e_code = zeros(1,code_n);
    for bm = 1:m_block
        for bn = 1:n_block
            for i = 1:8
                for j = 1:8
                    if (i+j==5)||(i+j==6)  %�е�Ƶ9��DCTϵ��
                        n_msg = n_msg + 1;
                        if n_msg<=code_n
                           yd = cover_DCT((bm-1)*8+i,(bn-1)*8+j);
                            tab = double(tab_m(i,j))/G;
                            dnum1 = round(yd/tab);
                            if mod(dnum1,2)==0
                                e_code(n_msg)=0;
                            else
                                e_code(n_msg)=1;
                            end
                        else
                            break;
                        end
                    end
                end
            end
            if n_msg>code_n break; end
        end
        if n_msg>code_n break; end
    end
    rand('seed',3); %��������
    SH = randperm(1*code_n);
    [e_code_Shuffle] = Image_Shuffle(e_code,SH);%����ͬ�������ҷ�������

%% У�������ȡ���������еľ���
    L=length(e_code_Shuffle);
    block=15;
    l_cover_embed=block*k*m;%C1����
    l_embed_check=3*k*m;
    % �ȴ�S3��ȡAdditional check code ����S2����
    stc_extract_Check_code_check = stc_extract(uint8(e_code_Shuffle(l_cover_embed+l_embed_check+1:L)'), l_check_check, 10);
    [decoded_check_bin] = rs_decode_klf(e_code_Shuffle(l_cover_embed+1:l_cover_embed+l_embed_check),stc_extract_Check_code_check,n,k,m);
    %�Ӿ�����S2����ȡ�����벢����S1
    stc_extract_Check_code = stc_extract(uint8(decoded_check_bin'), stc_n_msg_bits_check, 10);
     [decoded_msg_bin] = rs_decode_klf(e_code_Shuffle(1:l_cover_embed),stc_extract_Check_code,n,k,m);
    
%%  ��ȡ������Ϣ

    H = 10;
    stc_decoded_msg = stc_extract(uint8(decoded_msg_bin'), stc_n_msg_bits, H);
    stc_decoded_msg=stc_decoded_msg';

end
