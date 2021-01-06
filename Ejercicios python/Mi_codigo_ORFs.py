def complement(dna):
    complement = {
                "A":"T", 
                "C":"G", 
                "T":"A",
                "G":"C" 
	      }
    dna_complement = ''
    for i in range(0, len(dna)):
        dna_complement = dna_complement + complement[dna[i]]
        
    return dna_complement


def reverse_complement(dna):
    complement = {
                "A":"T", 
                "C":"G", 
                "T":"A",
                "G":"C" 
	      }
    dna_complement = ''
    for i in range(0, len(dna)):
        dna_complement = dna_complement + complement[dna[i]]
        
    dna_reverse_complement = dna_complement[::-1]
        
    return dna_reverse_complement

    
def transcription(dna):
    transcription = {
         	       "A":"A",
         	       "C":"C",
                   "T":"U",
                   "G":"G"
         	        }
    rna = ''
    for i in range(0, len(dna)):
        rna = rna + (transcription[dna[i]])
        
    return rna

def orfs(dna):
    aminoacidos = {
                  "UUU": "F",
                  "UUC": "F",
                  "UUA": "L",
                  "UUG": "L",
                  "UCU": "S",
                  "UCC": "S",
                  "UCA": "S",
                  "UCG": "S",
                  "UAU": "Y",
                  "UAC": "Y",
                  "UAA": "*",
                  "UAG": "*",
                  "UGU": "C",
                  "UGC": "C",
                  "UGA": "*",
                  "UGG": "W",
                  "CUU": "L",
                  "CUC": "L",
                  "CUA": "L",
                  "CUG": "L",
                  "CCU": "P",
                  "CCC": "P",
                  "CCA": "P",
                  "CCG": "P",
                  "CAU": "H",
                  "CAC": "H",
                  "CAA": "Q",
                  "CAG": "Q",
                  "CGU": "R",
                  "CGC": "R",
                  "CGA": "R",
                  "CGG": "R",
                  "AUU": "I",
                  "AUC": "I",
                  "AUA": "I",
                  "AUG": "M",
                  "ACU": "T",
                  "ACC": "T",
                  "ACA": "T",
                  "ACG": "T",
                  "AAU": "N",
                  "AAC": "N",
                  "AAA": "L",
                  "AAG": "L",
                  "AGU": "S",
                  "AGC": "S",
                  "AGA": "R",
                  "AGG": "R",
                  "GUU": "V",
                  "GUC": "V",
                  "GUA": "V",
                  "GUG": "V",
                  "GCU": "A",
                  "GCC": "A",
                  "GCA": "A",
                  "GCG": "A",
                  "GAU": "D",
                  "GAC": "D",
                  "GAA": "E",
                  "GAG": "E",
                  "GGU": "G",
                  "GGC": "G",
                  "GGA": "G",
                  "GGG": "G"
                }
    rna = transcription(dna)
    rev_rna = transcription(reverse_complement(dna))
    
    rfs = []
    
    rf1 = []
    rf2 = []
    rf3 = []
    rf4 = []
    rf5 = []
    rf6 = []
    
    rf1.append('')
    rf2.append('')
    rf3.append('')
    rf4.append('')
    rf5.append('')
    rf6.append('')
    
    rf1.append([])
    rf2.append([])
    rf3.append([])
    rf4.append([])
    rf5.append([])
    rf6.append([])
    
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        
        if len(codon) == 3:
            rf1[0] = rf1[0] + aminoacidos[rna[i:i+3]]
            rf1[1].append(i+1)
            rf1.append('positiva')
    rfs.append(rf1)
        
    for i in range(1, len(rna), 3):
        codon = rna[i:i+3]
        
        if len(codon) == 3:
            rf2[0] = rf2[0] + aminoacidos[rna[i:i+3]]
            rf2[1].append(i+1)
            rf2.append('positiva')
    rfs.append(rf2)
        
    for i in range(2, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) == 3:
            rf3[0] = rf3[0] + aminoacidos[rna[i:i+3]]
            rf3[1].append(i+1)
            rf3.append('positiva')
    rfs.append(rf3)
    
    for i in range(0, len(rev_rna), 3):
        codon = rev_rna[i:i+3]
        if len(codon) == 3:
            rf4[0] = rf4[0] + aminoacidos[rev_rna[i:i+3]]
            rf4[1].append(i+1)
            rf4.append('negativa')
    rfs.append(rf4)
    
    for i in range(1, len(rev_rna), 3):
        codon = rev_rna[i:i+3]
        if len(codon) == 3:
            rf5[0] = rf5[0] + aminoacidos[rev_rna[i:i+3]]
            rf5[1].append(i+1)
            rf5.append('negativa')
    rfs.append(rf5)
    
    for i in range(2, len(rev_rna), 3):
        codon = rev_rna[i:i+3]
        if len(codon) == 3:
            rf6[0] = rf6[0] + aminoacidos[rev_rna[i:i+3]]
            rf6[1].append(i+1)
            rf6.append('negativa')
    rfs.append(rf6)
            
    flag_orf = False
    orfs = []
    orf = ''    
    meta_rf = {}
    lista = []
    cont = 0

    for rf in rfs:
        
        for i in range(0, len(rf[0])):
            
            if rf[0][i] == 'M':
                flag_orf = True
                
            if flag_orf == True:
                orf = orf + rf[0][i]
                cont = cont + 1
                
                if cont == 1:
                    meta_rf['inicio'] = rf[1][i]
            
            if flag_orf == True and rf[0][i] == '*':
                flag_orf = False
                meta_rf['num_proteinas'] = cont
                meta_rf['num_bases'] = cont * 3
                meta_rf['paro'] = rf[1][i] + 2
                meta_rf['cadena'] = rf[2]
                
                if len(orf) > 100:
                    lista.append(orf)
                    lista.append(meta_rf)
                    orfs.append(lista)
                orf = '' 
                meta_rf = {}
                lista = []
                cont = 0

    
    num = len(orfs)
    
    print('ORFs:', num)
    
    for i in orfs:
        print(i[0][0:8]+'...'+str(i[0][-8:len(i[0])]),
              '['+str(i[1]['num_bases'])+' bp]',
              str(i[1]['inicio'])+'...'+str(i[1]['paro']),
              '-'+str(i[1]['cadena']))



    return orfs, num



dna = ""

f = open("covid19.txt", "r")
for x in f:
  dna = dna + x
dna = dna.replace("\n","")

ORFs = orfs(dna)


