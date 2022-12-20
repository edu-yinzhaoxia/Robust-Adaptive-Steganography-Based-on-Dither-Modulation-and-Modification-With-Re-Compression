function Test_Robustness_DMMR_JUNIWARD()
%% 基于压缩修改的测试代码
    clear all;
    clc;
%%  参数设置
    addpath(fullfile('jpeg_toolbox'));
    addpath(fullfile('STC'));
    cover_dir = ''; %载体jpeg图像所在文件夹
    recompress_dir='';if ~exist(recompress_dir,'dir'); mkdir(recompress_dir); end %中间图像所在的文件夹
    stego_dir = ''; if ~exist(stego_dir,'dir'); mkdir(stego_dir); end  %载密图像所在文件夹
    afterchannel_stego_dir = ''; if ~exist(afterchannel_stego_dir,'dir'); mkdir(afterchannel_stego_dir); end  %信道处理后载密图像所在文件夹     
    cover_num = 1000; %测试载体图像个数
    cover_QF = 65; %载体图像的质量因子
    attack_QF = 65; %模拟信道压缩的信道质量因子
    payload = 0.1; %嵌入率   
    bit_error_rate = zeros(1,cover_num); %记录测试图像的误码率
%%  消息嵌入    
    for i_img = 1:cover_num
        cover_Path = fullfile([cover_dir,'\',num2str(i_img),'.jpg']);  
        recompress_Path=fullfile([ recompress_dir,'\',num2str(i_img),'.jpg']);
        stego_Path = fullfile([stego_dir,'\',num2str(i_img),'.jpg']);    
        afterchannel_stego_Path = fullfile([afterchannel_stego_dir,'\',num2str(i_img),'.jpg']);   
        C_STRUCT = jpeg_read(cover_Path);
        C_COEFFS = C_STRUCT.coef_arrays{1};  
        C_QUANT = C_STRUCT.quant_tables{1}; %载体图像量化表
        nzAC = nnz(C_COEFFS) - nnz(C_COEFFS(1:8:end,1:8:end));
%  随机产生均匀分布的二进制原始秘密信息
        raw_msg_len = ceil(payload*nzAC);
        if raw_msg_len<=1 %对于BOSSbase数据集中存在不适合嵌入的不进行嵌入
            imwrite(imread(cover_Path),recompress_Path,'quality',cover_QF);
            imwrite(imread(recompress_Path),stego_Path,'quality',cover_QF);
            imwrite(imread(stego_Path),afterchannel_stego_Path,'quality',attack_QF);
            bit_error_rate(1,i_img)=0;
            continue;
        end
        raw_msg = round( rand(1,raw_msg_len) ); %原始秘密信息的行向量    
        n = 255; k = 251; m =8 ;    % RS编码参数 
        [rho1_P,rho1_M] = J_UNIWARD_D(cover_Path,1);%使用J_UNIWARD隐写算法计算DCT系数的嵌入代价
        [cover_lsb, change, rho] = DMMR(cover_Path, rho1_P, rho1_M, C_QUANT);%计算通过抖动调制计算载序列，修改距离，修改代价
        %  利用STC进行消息嵌入  
        [stc_n_msg_bits,stc_n_msg_bits_check,l_check_check] = stc_embedding(raw_msg, cover_Path, cover_lsb, rho, change,cover_QF, stego_Path,n,k,m,recompress_Path,C_QUANT,attack_QF);
%%  模拟信道压缩 
        imwrite(imread(stego_Path),afterchannel_stego_Path,'quality',attack_QF);    
      
%%  消息提取
        [stc_decoded_msg] = stc_extracting(afterchannel_stego_Path, stc_n_msg_bits,stc_n_msg_bits_check, C_QUANT,n,k,m,l_check_check);   

%%  计算每张图像的误码率      
        bit_error = double(raw_msg) - double(stc_decoded_msg);
        bit_error_number = sum(abs(bit_error));
        bit_error_rate(1,i_img) = bit_error_number/raw_msg_len;
        
%% 计算每张出错率
%  输出每张图像的误码率               
        fprintf('%s\n',['payload: ',num2str(payload),'  image_number: ',num2str(i_img),'  error_rate: ',num2str(bit_error_rate(1,i_img))]);  %
    end
%  输出所有图像的平均误码率
    ave_error_rate = mean(bit_error_rate);
    fprintf('%s\n',['payload: ',num2str(payload),'  ave_error_rate: ',num2str(ave_error_rate)]);  
  
end
