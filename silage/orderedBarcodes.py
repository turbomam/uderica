n = 20
from itertools import product
barcodes = [''.join(i) for i in product('ACGT', repeat = 8)]
for line in barcodes[:n]: print line
