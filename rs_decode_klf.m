function [decoded_msg_bin,Check_code_real] = rs_decode_klf(MSG1,MSG2,n,k,m)%n=255,k=253,m=8
L=length(MSG1);
L_dec=L/m;
n_block=L_dec/k;
func = @(x) (x(1)*2^(m-1) + x(2)*2^(m-2) + x(3)*2^(m-3) + x(4)*2^(m-4) + x(5)*2^(m-5)+x(6)*2^(m-6) + x(7)*2^(m-7) + x(8)*2^(m-8));

MSG1_dec=blkproc(MSG1,[1 m],func);%编码序列每8bit转换为10进制
Check_code_real = blkproc(MSG2',[1 m],func);  %将二进制校验码转换为10进制

Check_code_real_reshape=reshape(Check_code_real,[n-k,n_block]);
MSG1_reshape=reshape(MSG1_dec,[k,n_block]);
encoded_MSG=[MSG1_reshape;Check_code_real_reshape];
encoded_MSG=encoded_MSG';
encoded_MSG_gf=gf(encoded_MSG,m);
[decoded_MSG_gf,numerr,corrcode_gf]=rsdec(encoded_MSG_gf,n,k);%RS解码

decoded_msg= double(zeros(1,k*n_block));
corrcode=double(zeros(1,n*n_block));
encoded_MSG_1=double(zeros(1,n*n_block));

%% 将解码后的信息转换为数组
for i =1:n_block
    for j=1:k
        decoded_MSG_gf_x = double(decoded_MSG_gf.x);
        decoded_msg(1,(i-1)*k+j)=decoded_MSG_gf_x(i,j);
    end
end

%% 对于超过纠错能力的情况不予纠正
for i =1:n_block
    for j=1:n
        corrcode_gf_x = double(corrcode_gf.x);
        encoded_MSG_gf_x=double(encoded_MSG_gf.x);
        corrcode(1,(i-1)*n+j)=corrcode_gf_x(i,j);
        encoded_MSG_1(1,(i-1)*n+j)=encoded_MSG_gf_x(i,j);
    end
end

for i =1:n_block

    if numerr(i)==-1 %该片段超过纠错能力，不进行纠错
        decoded_msg(1,(i-1)*k+1:i*k)=MSG1_dec(1,(i-1)*k+1:i*k); 
    end
    if numerr(i)==(n-k)/2  %对于纠错数目等于纠错能力时可能存在的错误，采用每8bit出错个数不超过1个原则进行纠错，超过则不进行纠错
        b=abs(double(corrcode((i-1)*n+1:i*n))-double(encoded_MSG_1((i-1)*n+1:i*n)));
        c=find(b~=0);
        for ii=1:length(c)
            bin=dec2bin(b(c(ii)));
            d=sum(bin(:)=='1');
            if d>=2           %8bit出错个数大于1个，不进行纠错
                decoded_msg(1,(i-1)*k+1:i*k)=MSG1_dec(1,(i-1)*k+1:i*k); 
                
            end
        end
    end
end

%% 转换为二进制
decoded_msg_bin=zeros(1,k*n_block*m);
for i=1:k*n_block
    dec=decoded_msg(1,i);
    bin=dec2bin(dec,m);
    for j=1:m
        decoded_msg_bin(1,(i-1)*m+j)=str2num(bin(j));
    end
end

end







