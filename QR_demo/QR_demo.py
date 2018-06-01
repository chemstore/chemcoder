from collections import Counter
from chemcoder import ChemCoder, flags2binary, getFragmentLength

def readQR(filename):
    
    return ''.join(
        line.rstrip('\n').replace(' ', '0').replace('#', '1')
        for line in open(filename, 'r')
    )
    
# read the QR file as a bitstring
QR = readQR('QR_benzene.txt')
# print(QR)

flag_count = 15
polymer_length = 6

chemCoder = ChemCoder(flag_count=flag_count, polymer_length=polymer_length)
print(QR[:100])
fragments = chemCoder.encode(QR)
# print(len(fragments))
# print(fragments)
QR2 = chemCoder.decode(fragments)
# print(QR2)

# print(QR == QR2)

frequencyTable = Counter()
for fragment in fragments:
    for letter in fragment:
        frequencyTable[letter] += 1

# print list of fragments
for index, fragment in enumerate(fragments):
    # print('{:2d}: {:>6s}'.format(index + 1, fragment))
    print('{:2d} {} {:>6s}'.format(
        index + 1, 
        str(flags2binary(fragment, base=flag_count)).zfill(
            getFragmentLength(flag_count, polymer_length)
        ), 
        fragment
    ))
    
# print frequency table
print()
for letter, count in sorted(frequencyTable.items(), key=lambda x:(-x[1], x[0])):
    print('{}: {}'.format(letter, count))