function [stc_decoded_msg] = stc_extracting(afterchannel_stego_Path,stc_n_msg_bits,stc_n_msg_bits_check,tab_m,n,k,m,l_check_check)
%% 提取载密序列元素
    bits = 8;
    cover_spa = imread(afterchannel_stego_Path);
    cover_spa = double(cover_spa) - 2^(round(bits)-1);
    [xm,xn] = size(cover_spa);
    t = dctmtx(8);
    fun = @(xl) (t*xl*(t'));
    cover_DCT = blkproc(cover_spa,[8 8],fun);%分块DCT变换
    m_block = floor(xm/8);
    n_block = floor(xn/8);
    G = 1;
    n_msg = 0;
    code_n = m_block*n_block*9;  %中频9个DCT系数
    e_code = zeros(1,code_n);
    for bm = 1:m_block
        for bn = 1:n_block
            for i = 1:8
                for j = 1:8
                    if (i+j==5)||(i+j==6)  %中低频9个DCT系数
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
    rand('seed',3); %设置种子
    SH = randperm(1*code_n);
    [e_code_Shuffle] = Image_Shuffle(e_code,SH);%按照同样的置乱方法置乱

%% 校验码的提取和载密序列的纠错
    L=length(e_code_Shuffle);
    block=15;
    l_cover_embed=block*k*m;%C1长度
    l_embed_check=3*k*m;
    % 先从S3提取Additional check code 并对S2纠错
    stc_extract_Check_code_check = stc_extract(uint8(e_code_Shuffle(l_cover_embed+l_embed_check+1:L)'), l_check_check, 10);
    [decoded_check_bin] = rs_decode_klf(e_code_Shuffle(l_cover_embed+1:l_cover_embed+l_embed_check),stc_extract_Check_code_check,n,k,m);
    %从纠错后的S2中提取纠错码并纠错S1
    stc_extract_Check_code = stc_extract(uint8(decoded_check_bin'), stc_n_msg_bits_check, 10);
     [decoded_msg_bin] = rs_decode_klf(e_code_Shuffle(1:l_cover_embed),stc_extract_Check_code,n,k,m);
    
%%  提取秘密信息

    H = 10;
    stc_decoded_msg = stc_extract(uint8(decoded_msg_bin'), stc_n_msg_bits, H);
    stc_decoded_msg=stc_decoded_msg';

end
