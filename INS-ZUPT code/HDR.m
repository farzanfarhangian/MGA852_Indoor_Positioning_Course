function mk = HDR(saw_k,saw_k_1)
th=30*pi/180;
if abs(saw_k-saw_k_1)<=th
    mk=saw_k-saw_k_1;
else
    mk=0;
    %mk=saw_k;
end

end
