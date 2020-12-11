for f in /folder_with_fastq_files/*_1.fq.gz; do b=`basename $f`; c=(${b//_/ }); echo "${c[0]}_${c[1]}_${c[2]}_1i2"; kraken --thread $THREADS --fastq-input --gzip-compressed --preload -db $KRAKEN_DB_NAME --output /output_folder/"${c[0]}_${c[1]}_${c[2]}".report --paired /folder_with_fastq_files/"${c[0]}_${c[1]}_${c[2]}_1.fq.gz" /folder_with_fastq_files/"${c[0]}_${c[1]}_${c[2]}_2.fq.gz" --check-names; done

for f in /output_folder/*report; do echo `basename $f`; (kraken-report --db $KRAKEN_DB_NAME $f >$f.csv); done

