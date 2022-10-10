#import csv 


import pandas as pd
pd.set_option('display.max_columns', None)

here = '/Users/surasaksangkhathat/Desktop/VCFRIDER/'

def setpath(file, type):
    if type == 'snp':
        file = here+file+'.ann.snp.txt'
    elif type == 'indel':
        file = here+file+'.ann.indel.txt'
    return file

def vcfimport(path):
    vcf = pd.read_csv(path, delimiter="\t")
    vcf['GEN[*].GT'] = vcf['GEN[*].GT'].replace('1-Jan', '1/1')
    print(vcf.groupby(['ANN[*].EFFECT'])['ANN[*].EFFECT'].count())
    print('==================================================================')
    print(vcf.groupby(['ANN[*].IMPACT'])['ANN[*].IMPACT'].count())
    return vcf



def geneselect(vcf, genes = []):
	if len(genes) == 0:
		genes = list(map(str, input('Enter a list of factors to be analysed,separated by a comma: ').split(',')))
		genes = [i.strip() for i in genes]
		print(genes)
	else: 
		pass
	output = vcf.loc[vcf['ANN[*].GENE'].isin(genes)]
	print('Total number of variants captured',len(output))
	output_impact = output.groupby(['ANN[*].GENE','ANN[*].EFFECT'])['ANN[*].EFFECT'].count()
	print(output_impact)
	output = output[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GEN[*].GT',
       'GEN[*].AD[0]', 'GEN[*].AD[1]', 'GEN[*].DP', 'ANN[*].EFFECT', 'ANN[*].IMPACT',
       'ANN[*].CDNA_POS', 'ANN[*].AA_POS', 'ANN[*].GENE', 'ANN[*].BIOTYPE',
       'ANN[*].HGVS_C', 'ANN[*].HGVS_P', 'ANN[*].FEATUREID']]
	return output


def filtersnp(select):
    snplist = ['missense_variant', 'missense_variant&splice_region_variant', 'protein_protein_contact','splice_acceptor_variant&intron_variant','splice_donor_variant&intron_variant','start_lost',
	'start_lost&splice_region_variant','stop_gained','stop_lost','stop_lost&splice_region_variant', 'stop_retained_variant', 'start_lost&splice_region_variant', 'stop_gained&splice_region_variant', 'start_lost&splice_region_variant']
    output = select[(select['ANN[*].EFFECT'].isin(snplist))|(select['ANN[*].IMPACT'] == 'HIGH')]
    output = output.drop_duplicates(subset=['CHROM','POS','ANN[*].GENE'], keep='first')
    print('Total number of variants passed',len(output))
    output_impact = output.groupby(['ANN[*].GENE','ANN[*].EFFECT'])['ANN[*].EFFECT'].count()
    print(output_impact)
    return output
   
def filterindel(select):
    indellist = ['disruptive_inframe_deletion','frameshift_variant','frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant',
    'frameshift_variant&splice_region_variant', 'frameshift_variant&start_lost',
    'frameshift_variant&stop_gained','frameshift_variant&stop_lost',
    'start_lost&conservative_inframe_deletion','start_lost&disruptive_inframe_deletion',
    'stop_gained&disruptive_inframe_deletion','stop_gained&disruptive_inframe_insertion',
    'structural_interaction_variant', 'splice_acceptor_variant&intron_variant']
    output = select[(select['ANN[*].EFFECT'].isin(indellist))|(select['ANN[*].IMPACT']=='HIGH')]
    output = output.drop_duplicates(subset=['CHROM','POS','ANN[*].GENE'], keep='first')
    print('Total number of variants passed',len(output))
    output_impact = output.groupby(['ANN[*].GENE','ANN[*].EFFECT'])['ANN[*].EFFECT'].count()
    print(output_impact)
    return output


def export(vcf, name):
    vcf.to_csv(name+'.csv',index = False, sep=',', encoding='utf-8')

def prioritize(file, geneset):
    print(file)
    print('SNP filtering')
    path_snp = setpath(file, 'snp')
    df_snp = vcfimport(path_snp)
    df_snp_sel = geneselect(df_snp, geneset)
    print('FIltering ===============================================')
    df_snp_sel_filter = filtersnp(df_snp_sel)
    export(df_snp_sel_filter, file+'.snp.filtered')
    print ('========================================================')
    print('INDEL filtering')
    path_indel = setpath(file, 'indel')
    df_indel = vcfimport(path_indel)
    df_indel_sel = geneselect(df_indel, geneset)
    print('FIltering ===============================================')
    df_indel_sel_filter = filterindel(df_indel_sel)
    export(df_indel_sel_filter, file+'.indel.filtered')
     





#St_Jude gene sets for pediatric cancers

tyr_kinase = ['FLT3', 'ABL1', 'KIT', 'ALK', 'PDGFRA']
ras = ['NRAS', 'KRAS', 'PTPN11', 'NF1', 'BRAF', 'NF2', 'RIT1']
pi3k = ['PTEN', 'PIK3R1', 'AKI1', 'PIK3CD', 'FGFR1', 'PDFGRB', 'PIK3CA', 'TSC2', 'ATM', 'TSC1']
jak_stat = ['JAK3', 'JAK2', 'CRLF2', 'IL7R', 'CBL', 'PTPN2', 'STAT5B', 'JAK1', 'CSF3R', 'SH2B3', 'EPOR', 'XPO1']
notch = ['NOTCH', 'FBXW7'] 
wnt = ['CTNNB1', 'AMER1', 'APC']
mirna = ['DROSHA', 'DGCRB', 'XPO5']
myc = ['MYCN', 'MGA', 'MYC', 'MAX', 'TRRAP']
transcription = ['TAL1', 'IKZF1', 'PAX5', 'ETV6', 'RUNX1',
'WT1', 'LEF1', 'MYB', 'BCL11B', 'CBFB',
'TCF3', 'TLX3', 'LMO2', 'NUP98', 'EBF1',
'MLLT1', 'ERG', 'NKX2-1', 'TCF7', 'ZNF384',
'XBP1', 'TBL1XR1', 'GATA2', 'CNOT3', 'MED12',
'HOXA10', 'CEBPA', 'ELF1', 'GATA3', 'ZNF217',
'MEF2D', 'ZBTB7A', 'ZFP36L2', 'TAL2', 'NR3C1', 
'TSPYL2', 'EWSR1', 'GLIS2', 'LMO1', 'SIX1',
'NR3C2', 'UBTF', 'IKZF3', 'ELF4', 'IKZF2', 
'SIX2', 'GATA1', 'FLI1', 'DUX4', 'GATA4',
'GON4L', 'MEF2C', 'ZEB2', 'LYL1', 'FEV',
'ETS2'] 
cell_cycle = ['CDKN2A', 'TP53', 'RB1', 'CDKN1B', 'BTG1', 'TERT', 'CCND3', 'CDKN2B', 'NPM1', 'CDK6', 'CCND2', 'CDK4']
cohesin = ['STAG2', 'NIPBL', 'SMC3', 'PDS5B', 'STAG1', 'PDS5A', 'RAD21']
epigen = ['PHF6', 'KMT2A', 'CTCF', 'ATRX', 'CREBBP',
'KDM6A', 'EXH2', 'SETD2', 'SUZ12', 'ASXL2',
'EP300', 'TOX', 'KMT2D', 'SMARCA4', 'BCOR',
'ASXL1', 'EED', 'BCORL1', 'ARID2', 'BAZ1A',
'TET2', 'ARID1A', 'WHSC1', 'ATF7IP', 'MLLT4',
'NCOR1', 'KDM5A', 'HDAC7', 'NSD1', 'KMT2E',
'KMT2C', 'CHD4', 'DNMT3A', 'KDM5C', 'HIST1H3C',
'PBRM1', 'CHD7', 'TET3', 'KMT2B', 'HDA39']
st_jude = tyr_kinase + ras + pi3k + jak_stat + notch + wnt + mirna + myc +transcription + cell_cycle + cohesin + epigen

drugable = ['AKT1', 'ALK', 'ARAF', 'ARID1A', 'ATM', 'BARD', 'BRAF'
'BRCA1', 'BRCA2', 'BRIP', 'CDK4', 'CDK12', 'CDKN2A', 'CHEK1', 'CHEK2'
'EGFR', 'ERBB2', 'ERCC2', 'FGFR3', 'FLI1', 'HRAS'
'IDH1', 'KDM6A', 'KIT', 'KRAS', 'MAP2K1', 'MDM2',
'MET', 'MTOR', 'NF1', 'NRAS', 'NRG1', 'NTRK1', 'NTRK2', 'NTRK3'
'PALB2', 'PDGFB', 'PDGFRA', 'PIK3A', 'PTCH1', 'PTEN', 'RAD51B'
'RAD51C', 'RAD51D', 'RAD54L', 'RET', 'ROS1', 'SMARCB1', 'STK11', 'TSC1', 'TSC2']

tum = ['A008', 'A012', 'A015', 'A018', 'A020', 'A021', 'A025', 'A028', 'A031', 'A035', 'A040', 'A044']  
nb_list = ['A015', 'A018', 'A020', 'A021', 'A025', 'A028', 'A031', 'A035', 'A040', 'A044', 'A045']


def batchfilter(tum, geneset):
    for i in range(len(tum)):
        prioritize(tum[i], geneset)
        
def genedict(filelist, folder):
    gene_dict = {}
    for i in filelist:
        read_snp = pd.read_csv(here+folder+'/'+i+'.snp.filtered.csv')
        gene_list_i_snp = read_snp['ANN[*].GENE'].to_list()
        read_indel = pd.read_csv(here+folder+'/'+i+'.indel.filtered.csv')
        gene_list_i_indel = read_indel['ANN[*].GENE'].to_list()
        gene_list_i_mix = gene_list_i_snp + gene_list_i_indel
        for entry in gene_list_i_mix:
            if entry in (gene_dict):
                gene_dict[entry] += 1
            else:
                    gene_dict[entry] = 1
    gene_dict = {str(key): val for key, val in gene_dict.items()}
    for key in sorted(gene_dict.keys()):
        print(key , " :: " , gene_dict[key])
        print('')
    print(sorted(gene_dict.items(), key=lambda x: x[1], reverse = True))
    return gene_dict

def genedict_unique(filelist, folder):
    gene_dict_case = {}
    for i in filelist:
        read_snp = pd.read_csv(here+folder+'/'+i+'.snp.filtered.csv')
        gene_list_i_snp = read_snp['ANN[*].GENE'].to_list()
        read_indel = pd.read_csv(here+folder+'/'+i+'.indel.filtered.csv')
        gene_list_i_indel = read_indel['ANN[*].GENE'].to_list()
        gene_list_i_indel = gene_list_i_indel + gene_list_i_snp
        gene_list_i_indel = set(gene_list_i_indel)
        for entry in gene_list_i_indel:
            if entry in (gene_dict_case):
                gene_dict_case[entry] += 1
            else:
                    gene_dict_case[entry] = 1
    gene_dict_case = {str(key): val for key, val in gene_dict_case.items()}
    for key in sorted(gene_dict_case.keys()):
        print(key , " :: " , gene_dict_case[key], 'cases ,percentage of variants = ', gene_dict_case[key]*100/len(filelist),'%')
        print('')
    print(sorted(gene_dict_case.items(), key=lambda x: x[1], reverse = True))
    return gene_dict_case

def matrix(geneset, folder, file):
    df = pd.DataFrame({'Genelist':geneset})
    A_list = []
    read_snp = pd.read_csv(here+folder+'/'+file+'.snp.filtered.csv')
    gene_list_i_snp = read_snp['ANN[*].GENE'].to_list()
    read_indel = pd.read_csv(here+folder+'/'+file+'.indel.filtered.csv')
    gene_list_i_indel = read_indel['ANN[*].GENE'].to_list()
    for i in geneset:
        if i in gene_list_i_snp or i in gene_list_i_indel:
            A_list.append('1')
        else:
            A_list.append('0')
    namae = str(file)
    df[namae] = A_list
    return df 

def variantmatrix(geneset, folder, filelist):
    df = pd.DataFrame({'Genelist':geneset})
    for i in filelist:
        df2 = matrix(geneset, folder, i)
        df = pd.concat([df,df2[str(i)]], axis = 1)
    return df
                