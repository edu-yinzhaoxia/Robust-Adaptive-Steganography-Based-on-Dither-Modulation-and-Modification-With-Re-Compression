function [l,l_check,l_check_check] = stc_embedding(msg,coverPath,cover_lsb,rho,change,cover_QF,stegoPath,n,k,m,recompress_Path,tab_m,attack_QF)%用于修改参数
%% 设置信息嵌入和校验码嵌入的分界点
L=length(cover_lsb);
block=15;
l_cover_embed=block*k*m;%C1的长度
l=length(msg);

%% 将载体序列置乱
rand('seed',3); %设置种子
SH = randperm(1*L);
[cover_lsb_Shuffle] = Image_Shuffle(cover_lsb,SH);
[rho_Shuffle]=Image_Shuffle(rho,SH);

%% STCs 编码
H = 10;
[min_cost, stc_msg] = stc_embed(uint8(cover_lsb_Shuffle(1:l_cover_embed)'),  uint8(msg'),rho_Shuffle(1:l_cover_embed), H); % 嵌入信息    
stc_extract_msg2 = stc_extract(stc_msg, l, H); % 提取信息

%% 验证STC解码是否正常工作
if all(uint8(msg) == stc_extract_msg2')
    disp('Message can be extracted by STC3 correctly.');
else
    error('Some error occured in the extraction process of STC3.');
end

%% RS编码
l_embed_check=3*k*m;%C2的长度
[Check_code]=rs_encode_klf(stc_msg,n,k,m);%RS编码
l_check=length(Check_code);
[min_cost, stc_msg_check] = stc_embed(uint8(cover_lsb_Shuffle(l_cover_embed+1:l_cover_embed+l_embed_check)'),  uint8(Check_code'),rho_Shuffle(l_cover_embed+1:l_cover_embed+l_embed_check), H); % embed message    量化索引调制机制
stc_extract_check = stc_extract(stc_msg_check, l_check, H);
%% 验证STC解码是否正常工作
if all(uint8(Check_code) == stc_extract_check')
    disp('Message can be extracted by STC3 correctly.');
else
    error('Some error occured in the extraction process of STC3.');
end

%% RS编码
[Additional_check_code ]=rs_encode_klf(stc_msg_check,n,k,m);%生成Additional check code 
l_check_check=length(Additional_check_code);
[min_cost, stc_msg_check_check] = stc_embed(uint8(cover_lsb_Shuffle(l_cover_embed+l_embed_check+1:L)'),  uint8(Additional_check_code'),rho_Shuffle(l_cover_embed+l_embed_check+1:L), H); % embed message    量化索引调制机制
stc_extract_check_check = stc_extract(stc_msg_check_check, l_check_check, H);
%% 验证STC解码是否正常工作
if all(uint8(Additional_check_code) == stc_extract_check_check')
    disp('Message can be extracted by STC3 correctly.');
else
    error('Some error occured in the extraction process of STC3.');  
end
%% 载密序列的整合
cover_lsb_Shuffle(1:l_cover_embed)=stc_msg;
cover_lsb_Shuffle(l_cover_embed+1:l_cover_embed+l_embed_check)=stc_msg_check;
cover_lsb_Shuffle(l_cover_embed+l_embed_check+1:L)=stc_msg_check_check;
[cover_lsb_change] = Image_ReShuffle(cover_lsb_Shuffle,SH);

%%  DCT变换
code = cover_lsb_change;
code_n = length(cover_lsb);
bits = 8;     
cover_spa = imread(coverPath);
cover_spa = double(cover_spa) - 2^(round(bits)-1);
[xm,xn] = size(cover_spa);
t = dctmtx(8);  %产生DCT矩阵
fun = @(xl) (t*xl*(t'));
cover_DCT = blkproc(cover_spa,[8 8],fun);%分块DCT变换
m_block = floor(xm/8);
n_block = floor(xn/8);

%% 嵌入过程
G = 1;
n_msg = 0;
for bm = 1:m_block
    for bn = 1:n_block
        for i = 1:8
            for j = 1:8
                if (i+j==5)||(i+j==6)  %中频21个DCT系数(i+j==6)||(i+j==7)%
                    n_msg = n_msg + 1;
                    if n_msg<=code_n
                        yd = cover_DCT((bm-1)*8+i,(bn-1)*8+j); %非量化的DCT系数
                        if code(n_msg) ~= cover_lsb(n_msg)   % 实际修改操作
                            yd = yd + change(n_msg); 
                        else
                            yd = cover_DCT((bm-1)*8+i,(bn-1)*8+j); %对于code(n_msg) == cover_lsb(n_msg)，即非嵌入位置保持不变
                        end
                        cover_DCT((bm-1)*8+i,(bn-1)*8+j) = yd; 
                    else
                        break;
                    end

                end
            end
        end
        if n_msg>code_n
            break; end
    end
    if n_msg>code_n
        break; end
end


%% 修改阶段 修改不稳定的DCT系数
for i=1:2
    cover_spa = blkproc(cover_DCT,[8 8],'P1*x*P2',t',t);
    cover_spa = cover_spa + double(2^(bits-1));
    cover_spa = uint8(cover_spa);
    imwrite(cover_spa,recompress_Path,'quality',attack_QF);
    bits = 8;     
    cover_spa = imread(recompress_Path);
    cover_spa = double(cover_spa) - 2^(round(bits)-1);
    t = dctmtx(8);  %产生DCT矩阵
    fun = @(xl) (t*xl*(t'));
    cover_DCT = blkproc(cover_spa,[8 8],fun);%分块DCT变换
    n_msg = 0;
    code_n = m_block*n_block*9;  %中低频9个DCT系数
    e_code = zeros(1,code_n);
    for bm = 1:m_block
        for bn = 1:n_block
            for i = 1:8
                for j = 1:8
                    if (i+j==5)||(i+j==6)  %中频9个DCT系数 提取载密序列
                        n_msg = n_msg + 1;
                         if n_msg<=code_n
                            yd = cover_DCT((bm-1)*8+i,(bn-1)*8+j);
                            tab_q = double(tab_m(i,j))/G;
                            dnum1 = round(yd/tab_q);
                            if mod(dnum1,2)==0
                                e_code(n_msg)=0;
                            else
                                e_code(n_msg)=1;
                            end
                        else
                            break;
                        end
                        if e_code(n_msg)~=cover_lsb_change(n_msg)%将提取的载密序列和原载密序列对比，如果出错就进行对应非量化DCT系数的修改
                            dnum1 = round(yd/tab_q);
                             if mod(dnum1,2)==0
                                dnum2 = floor(yd/tab_q);
                                if mod(dnum2,2)==1
                                    yd = dnum2*tab_q;%减去这么多
                                else
                                    yd = (dnum2+1)*tab_q;%加上这么多
                                end                        
                             else
                                dnum2 = floor(yd/tab_q);
                                if mod(dnum2,2)==1
                                    yd = (dnum2+1)*tab_q;%加上这么多
                                else
                                    yd = dnum2*tab_q;%加上这么多
                                end
                             end  
                            cover_DCT((bm-1)*8+i,(bn-1)*8+j)=yd;
                        end
                    end
                end
            end
            if n_msg>code_n break; end
        end
        if n_msg>code_n break; end
    end  
end
%% 生成载密图像
cover_spa = blkproc(cover_DCT,[8 8],'P1*x*P2',t',t);
cover_spa = cover_spa + double(2^(bits-1));
cover_spa = uint8(cover_spa);
imwrite(cover_spa,stegoPath,'quality',cover_QF);

end