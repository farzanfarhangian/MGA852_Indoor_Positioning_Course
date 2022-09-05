function yaw_new=corridorfunc(yaw_v,zupt)
k=length(yaw_v);
w=390;
x=sort(find(zupt==1),'descend');
zupt_k=x(1);  zupt_k_1=x(2);
yaw_new = yaw_v(k);
if zupt(k)==1
    if abs(yaw_v(k)-yaw_v(zupt_k))<=w
        yaw_new = yaw_v(k-1);
    end
end
        
end