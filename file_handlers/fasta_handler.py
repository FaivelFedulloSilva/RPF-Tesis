from libs.DNAString import DNAString
from libs.DNAStringSet import DNAStringSet


class FastaHandler:
    def __init__(self, path) -> None:
        def read_fasta(path: str):
            referenceDNA = DNAStringSet()
            current_reference = []
            current_id = None
            with open(path, 'r') as file:
                line = file.readline()
                while line != '':
                    if line[0] == '>':
                        if current_id != None:
                            referenceDNA.add_sequence(current_id, DNAString(' '.join(current_reference)))
                        current_id = line[1:].split(' ')[0]
                        current_reference= []
                    else:
                        current_reference.append(line.strip())
                    line = file.readline()

        self.reference = read_fasta(path)

    def get_reference(self):
        return self.reference 


      

# faHandler = FastaHandler(r'./Data/reference/hg38.fa')