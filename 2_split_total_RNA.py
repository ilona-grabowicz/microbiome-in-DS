f=open('./files_combined_SK_total_paired.csv', 'r').readlines()
for level in ['P', 'C', 'O', 'F', 'G', 'S', '-']:	
	h=open('split_file_combined_total_RNA_'+level+'.csv', 'w')
	h.write(f[0])
	for line in f[1:]:
		line=line.strip().split(',')
		print (line)
		if line[1]==level:
			h.write(','.join(line)+'\n')
