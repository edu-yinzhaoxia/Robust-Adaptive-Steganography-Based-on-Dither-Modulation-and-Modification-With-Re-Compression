function [Check_code] = rs_encode_klf(MSG,n,k,m)%n=255,k=253,m=8
%载体流的切分
L=length(MSG);
L_dec=L/m;
n_block=L_dec/k;
func = @(x) (x(1)*2^(m-1) + x(2)*2^(m-2) + x(3)*2^(m-3) + x(4)*2^(m-4) + x(5)*2^(m-5)+x(6)*2^(m-6) + x(7)*2^(m-7) + x(8)*2^(m-8));

MSG_dec = blkproc(MSG',[1 m],func);%编码序列每8bit转换为10进制
MSG_reshape=reshape(MSG_dec,[k,n_block]);%将得到的序列按照限制长度切分
MSG_reshape=MSG_reshape';
MSG_reshape_gf=gf(MSG_reshape,m);
encoded_MSG_reshape_gf=rsenc(MSG_reshape_gf,n,k);%RS编码

Check_code_len=(n-k)*n_block*m;
Check_code = zeros(1,Check_code_len);
index=1;
%将得到的十进制校验码转换为2进制序列
for i=1:n_block
    for j=k+1:n
        encoded_MSG_reshape_gf_x=double(encoded_MSG_reshape_gf.x);
        Check_code_each=dec2bin(encoded_MSG_reshape_gf_x(i,j),m);       
        for kk=1:m
            Check_code(1,index+kk-1) = double(str2num(Check_code_each(kk)));   
        end
        index=index+m;
        
    end
end

end
       
        
        
        
        
        
        


