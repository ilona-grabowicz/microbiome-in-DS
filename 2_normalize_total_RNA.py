f=open('./files_combined_SK_total_paired.csv', 'r').readlines()
for level in ['P', 'C', 'O', 'F', 'G', 'S', '-']:	
	lista_max=[]
	for i in range(len(f[0].strip().split(',')[2:-1])):
		lista_max.append(0)

	for line in f[3:]:
		line=line.strip().split(',')
		#print(line)
		if line[1]==level:
			for i in range(len(line[2:-1])):
				lista_max[i]+=int(line[2+i])	
	maximum=max(lista_max)

	h=open('split_file_combined_normalized_total_RNA_'+level+'.csv', 'w')
	h.write(f[0])
	for line in f[1:]:
		line=line.strip().split(',')
		if line[1]==level:
			h.write('%s,%s,' % (line[0], line[1]))
			roboczy=[]
			for i in range(len(line[2:-1])):
				print(line[2+i])
				roboczy.append(str(round(int(line[2+i])*maximum/lista_max[i])))
			h.write(','.join(roboczy)+'\n')	

	#print(lista_max)
	#print(maximum)
