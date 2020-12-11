import os

def prog():
	
	p = open('files_combined_SK_Biagi.csv', 'w')		
	p.write('taxon,level,')
	dict_taxon = {}
	all_taxa = []
	
#MAKING LIST OF ALL TAXA	
	for plik in os.listdir('./seed_kraken_results_Biagi_et_al/'):
		f=open('./seed_kraken_results_Biagi_et_al/'+plik).readlines()
		for line in f:
			line=line.strip().split('\t')
			if [line[5].strip(), line[3]] not in all_taxa:
				all_taxa.append([line[5].strip(), line[3]])
		p.write('%s,' % (plik[:-15]))
		print(plik[:-15])	
	p.write('\n')
	
#MAKING DICT WITH ALL FILES AND READ NUMBERS PER TAXON					
	for plik in os.listdir('./seed_kraken_results_Biagi_et_al/'):
		f=open('./seed_kraken_results_Biagi_et_al/'+plik).readlines()
		for line in f:
			line=line.strip().split('\t')	
						
			dict_taxon[(plik[:-15], line[5].strip())] = int(line[1])
	print (dict_taxon)
#WRITING DATA TO A FILE			
	for taxon in all_taxa:
		p.write('%s,%s,' % (taxon[0], taxon[1]))
		for plik in os.listdir('./seed_kraken_results_Biagi_et_al/'):
			if (plik[:-15], taxon[0]) in dict_taxon.keys():
				p.write('%i,' % (dict_taxon[(plik[:-15], taxon[0])]))
			else:
				
				p.write('0,')
		p.write('\n')
				
prog()


