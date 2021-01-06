# Código personal para leer archivos GenBank, encontrar su genoma y sus Genes (DNA, RNA, Proteinas)
# Utilicé una combinación de funciones y clases


def buscar(elemento,cadena):
    pos = cadena.find(elemento)
    tam = len(elemento)
    dato = ''
    for i in range((pos + tam + 2), len(cadena)):
        if cadena[i] == '\"':
            return dato
        dato = dato + cadena[i]

    return dato

    
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


def traslation(rna):
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
    protein = ''

    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) == 3:
            protein = protein + aminoacidos[rna[i:i+3]]
    return protein

    
class Genoma():   
    def __init__(self,file):
        self.file = file
        data = ''
        # info = {}
        self.genes = []
        f = open(self.file, 'r')
        for i in f:
          data = data + i
          
        # Extraer la Secuencia-------------------------------------
        self.dna = ''
        seq = ''
        for i in range(data.find('ORIGIN'), len(data)):
            seq = seq + data[i]
        
        for i in range(len(seq)):
            if seq[i] == 'a':
                self.dna = self.dna + 'A'
            if seq[i] == 't':
                self.dna = self.dna + 'T'
            if seq[i] == 'c':
                self.dna = self.dna + 'C'
            if seq[i] == 'g':
                self.dna = self.dna + 'G'
        
        # Extraer Información del genoma----------------------------
        features = ''
        if data.find('3\'UTR'):
            for i in range(data.find('FEATURES '), data.find('3\'UTR')):
                features = features + data[i]
        else:
            for i in range(data.find('FEATURES '), data.find('ORIGIN')):
                features = features + data[i]
    
        pos_genes = []
        posicion = 0
        text_genes = []
        while posicion != -1:
            posicion = features.find(' gene ',posicion)
            if posicion != -1:
                pos_genes.append(posicion)
                posicion += 1
        
        for i in range(len(pos_genes)):
            if i < (len(pos_genes)-1):
                text_genes.append(features[pos_genes[i]:pos_genes[i+1]])
            else:
                text_genes.append(features[pos_genes[i]:len(features)])
    
        for i in text_genes:
            dict = {}
            i = i.replace(' ','')
            i = i.replace('\n','')
    
            dict['gen'] = buscar('/gene',i)
            dict['locus'] = buscar('/locus_tag',i)
            dict['producto'] = buscar('/product',i)
            dict['id_protein'] = buscar('/protein_id',i)
            dict['inicio'] = str(i[(i.find('gene')+4):i.find('.')])
            dict['paro'] = str(i[(i.find('.')+2):i.find('/')])
            dict['seq'] = self.dna[int(dict['inicio'])-1:int(dict['paro'])]
    
            self.genes.append(dict)
        
        print('ORGANISMO: ' + buscar('/organism',features))
        print('MOLÉCULA: ' + buscar('/mol_type',features))
        print('TAMAÑO :', len(self.dna), 'bp')
        print('AISLADA EN : ' + buscar('/isolate',features))
        print('COLECTADA: ' + buscar('/collection_date',features))
    

    
    def seq(self):
        print(self.dna)
    
    
    
    def g(self,gen):
        for i in range(len(self.genes)) :
            if gen == self.genes[i]['gen']:
                return self.genes[i]
                

class Gen():
    def __init__(self,dict):
        self.inicio = dict['inicio']
        self.paro = dict['paro']
        self.seq = dict['seq']
        self.gen = dict['gen']
        self.producto = dict['producto']
        self.seq_rna = transcription(self.seq)
        self.protein = traslation((self.seq_rna))
        
        print(self.seq[0:8]+'...'+self.seq[-8:len(self.seq)],
                  '['+str(int(self.paro)-int(self.inicio))+' bp]',
                  str(self.inicio)+'...'+str(self.paro),
                  ' ('+str(self.gen)+') ',
                  '--> '+str(self.producto))
        

        
print('\n')
print('------------------------')

genoma = Genoma('sars-cov2.gb')
print('\n')
print('------------------------')

print('Datos del Gen de la proteina E:')
E = Gen(genoma.g('E'))
print('\n')
print('------------------------')

print('Secuencia del DNA de la proteina E:')
DNA = E.seq
print(DNA)
print('\n')
print('------------------------')

print('Secuencia del RNA de la proteina E:')
RNA = transcription(E.seq)
print(RNA)
print('\n')
print('------------------------')

print('Secuencia de la proteina E:')
Protein = traslation(RNA)
print(Protein)