function image_linear = ImageLinearFunction( image,k)
%对灰度值进行线性拉伸
%image: 被拉伸的影像
%k: 拉伸系数 一般为0.005

[m,n] = size(image);

%统计灰度直方图
bins = 0:intmax('uint16');
hist_count = histc(image(:),bins);

%清除临时变量
clear bins;
%直方图裁剪，计算左值和右值
[l_val,r_val] = calcu_val(hist_count,m*n,k);
%根据左值和右值对数据进行转换
image_linear = (image-l_val).*(255/(r_val-l_val));
image_linear = uint8(round(image_linear));

end

function [l_val,r_val] = calcu_val(hist_count,image_size,k)
%计算线性拉伸时对应的左值和右值

pix_cut = image_size*k;

%计算左值
temp = 0;
for i = 1:length(hist_count)
    temp = temp + hist_count(i);
    if temp >= pix_cut
        l_val = i;
        break;
    end
end
%计算右值
temp = 0;
for i = length(hist_count):-1:1
   temp = temp+ hist_count(i);
   if temp >= pix_cut
       r_val = i;
       break;
   end
end


end