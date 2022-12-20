function Test_Robustness_DMMR_JUNIWARD()
%% ����ѹ���޸ĵĲ��Դ���
    clear all;
    clc;
%%  ��������
    addpath(fullfile('jpeg_toolbox'));
    addpath(fullfile('STC'));
    cover_dir = ''; %����jpegͼ�������ļ���
    recompress_dir='';if ~exist(recompress_dir,'dir'); mkdir(recompress_dir); end %�м�ͼ�����ڵ��ļ���
    stego_dir = ''; if ~exist(stego_dir,'dir'); mkdir(stego_dir); end  %����ͼ�������ļ���
    afterchannel_stego_dir = ''; if ~exist(afterchannel_stego_dir,'dir'); mkdir(afterchannel_stego_dir); end  %�ŵ����������ͼ�������ļ���     
    cover_num = 1000; %��������ͼ�����
    cover_QF = 65; %����ͼ�����������
    attack_QF = 65; %ģ���ŵ�ѹ�����ŵ���������
    payload = 0.1; %Ƕ����   
    bit_error_rate = zeros(1,cover_num); %��¼����ͼ���������
%%  ��ϢǶ��    
    for i_img = 1:cover_num
        cover_Path = fullfile([cover_dir,'\',num2str(i_img),'.jpg']);  
        recompress_Path=fullfile([ recompress_dir,'\',num2str(i_img),'.jpg']);
        stego_Path = fullfile([stego_dir,'\',num2str(i_img),'.jpg']);    
        afterchannel_stego_Path = fullfile([afterchannel_stego_dir,'\',num2str(i_img),'.jpg']);   
        C_STRUCT = jpeg_read(cover_Path);
        C_COEFFS = C_STRUCT.coef_arrays{1};  
        C_QUANT = C_STRUCT.quant_tables{1}; %����ͼ��������
        nzAC = nnz(C_COEFFS) - nnz(C_COEFFS(1:8:end,1:8:end));
%  ����������ȷֲ��Ķ�����ԭʼ������Ϣ
        raw_msg_len = ceil(payload*nzAC);
        if raw_msg_len<=1 %����BOSSbase���ݼ��д��ڲ��ʺ�Ƕ��Ĳ�����Ƕ��
            imwrite(imread(cover_Path),recompress_Path,'quality',cover_QF);
            imwrite(imread(recompress_Path),stego_Path,'quality',cover_QF);
            imwrite(imread(stego_Path),afterchannel_stego_Path,'quality',attack_QF);
            bit_error_rate(1,i_img)=0;
            continue;
        end
        raw_msg = round( rand(1,raw_msg_len) ); %ԭʼ������Ϣ��������    
        n = 255; k = 251; m =8 ;    % RS������� 
        [rho1_P,rho1_M] = J_UNIWARD_D(cover_Path,1);%ʹ��J_UNIWARD��д�㷨����DCTϵ����Ƕ�����
        [cover_lsb, change, rho] = DMMR(cover_Path, rho1_P, rho1_M, C_QUANT);%����ͨ���������Ƽ��������У��޸ľ��룬�޸Ĵ���
        %  ����STC������ϢǶ��  
        [stc_n_msg_bits,stc_n_msg_bits_check,l_check_check] = stc_embedding(raw_msg, cover_Path, cover_lsb, rho, change,cover_QF, stego_Path,n,k,m,recompress_Path,C_QUANT,attack_QF);
%%  ģ���ŵ�ѹ�� 
        imwrite(imread(stego_Path),afterchannel_stego_Path,'quality',attack_QF);    
      
%%  ��Ϣ��ȡ
        [stc_decoded_msg] = stc_extracting(afterchannel_stego_Path, stc_n_msg_bits,stc_n_msg_bits_check, C_QUANT,n,k,m,l_check_check);   

%%  ����ÿ��ͼ���������      
        bit_error = double(raw_msg) - double(stc_decoded_msg);
        bit_error_number = sum(abs(bit_error));
        bit_error_rate(1,i_img) = bit_error_number/raw_msg_len;
        
%% ����ÿ�ų�����
%  ���ÿ��ͼ���������               
        fprintf('%s\n',['payload: ',num2str(payload),'  image_number: ',num2str(i_img),'  error_rate: ',num2str(bit_error_rate(1,i_img))]);  %
    end
%  �������ͼ���ƽ��������
    ave_error_rate = mean(bit_error_rate);
    fprintf('%s\n',['payload: ',num2str(payload),'  ave_error_rate: ',num2str(ave_error_rate)]);  
  
end
