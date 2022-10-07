BEGIN{FS=OFS="\t"}
{count[$4] += 1}
END{for(i=1;i<2000;i++) { count[i] = count[i] > 0 ? count[i] : 0; print i,count[i] } }
