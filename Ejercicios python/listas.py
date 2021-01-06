print('aqui trabajaremos con listas')

adnseq  = ["A", "C", "T", "G", "A", "T", "G", "T", "A", "C"]
adnseq2 = ["T", "G", "A", "T", "G"]

n = len(adnseq)
print (adnseq)
print (" tamano = ",n)
print (" primera posicion  = ",adnseq[0])
print (" ultima posicion= ",adnseq[n-1])
print (" primeras tres posiciones = ",adnseq[0:3])
print (" ultimas 5 posiciones = ",adnseq[-5:])

print ("a qui inseertamos un elemento...")
adnseq.append("T")
n = len(adnseq)
print (adnseq)
print (" tamano = ",n)

print ("remplaza la primera posicion en T")
adnseq[0] = "T"
print (adnseq)

print ("borrar la ultima posicion")
adnseq.pop(n-1)
print (adnseq)

print (" nueva secuancia = seq1 + seq2")

newadnseq = adnseq + adnseq2

print (adnseq)
print (adnseq2)
print (newadnseq)
print (" tamano = ", len(newadnseq))