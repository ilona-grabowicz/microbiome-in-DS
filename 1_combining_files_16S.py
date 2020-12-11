import os

def prog():
	
	p = open('files_combined_SK_16S.csv', 'w')		
	p.write('taxon,level,')
	dict_taxon = {}
	all_taxa = []
	
# MAKING DICT WITH ALL FILES AND READ NUMBERS PER TAXON	and MAKING LIST OF ALL TAXA	
	for plik in os.listdir('./seed_kraken_results_16S/'):
		f=open('./seed_kraken_results_16S/'+plik).readlines()
		for line in f:
			line=line.strip().split('\t')
			dict_taxon[(plik[:-20], line[5].strip())] = int(line[1])
			if [line[5].strip(), line[3]] not in all_taxa:
				all_taxa.append([line[5].strip(), line[3]])
		p.write('%s,' % (plik[:-20]))
	p.write('\n')
		
# WRITING DATA TO A FILE			
	for taxon in all_taxa:
		p.write('%s,%s,' % (taxon[0], taxon[1]))
		for plik in os.listdir('./seed_kraken_results_16S/'):
			if (plik[:-20], taxon[0]) in dict_taxon.keys():
				p.write('%i,' % (dict_taxon[(plik[:-20], taxon[0])]))
			else:
				
				p.write('0,')
		p.write('\n')
				
prog()


