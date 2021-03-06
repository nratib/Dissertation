###compile gd files into one folder###
for i in /Users/ratib/Desktop/breseq_data/*.fq.gz
do
	cp $i/output/output.gd /Users/ratib/Desktop/breseq_data/gd/${i##*/}.gd
done

###print all file paths to paste into gdtools###
for i in /Users/ratib/Desktop/breseq_T3clones/gd/*.gd
do
	echo $i '\'
done

gdtools COMPARE -o /Users/ratib/Desktop/breseq_T1clones/compare.tsv -f TSV -r /Users/ratib/zk126/ZK1142.gb /Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_20.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_21.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_22.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_23.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_C_24.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_20.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_21.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_22.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_23.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0010_S3_O_24.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_20.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_21.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_22.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_23.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_C_24.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_20.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_21.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_22.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_23.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0020_L3_M_24.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_20.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_21.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_22.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_23.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_24.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_25.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_26.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_27.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_28.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_29.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_30.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_31.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_32.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_33.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_35.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_36.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_37.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_38.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_39.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_40.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_41.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_42.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_43.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_44.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_45.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_46.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_47.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0030_L3_O_48.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_C_20.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_M_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_M_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_M_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_M_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_M_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_M_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_M_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_M_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/0060_L3_O_20.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_C_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_12.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_20.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_21.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_22.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_23.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_24.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_25.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_26.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_27.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_28.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_29.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_30.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_31.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/1050_L3_M_32.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_C_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_C_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_C_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_C_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_C_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_C_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_C_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_C_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_01.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_02.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_03.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_04.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_05.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_06.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_07.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_08.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_09.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_10.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_11.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_13.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_14.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_15.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_16.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_17.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_18.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_19.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_21.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_22.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_23.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_24.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_25.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_26.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_27.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_28.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_29.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_30.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_31.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_32.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_33.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_34.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_35.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_36.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_37.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_38.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_39.fq.gz.gd \
/Users/ratib/Desktop/breseq_T3clones/gd/3800_L3_M_40.fq.gz.gd