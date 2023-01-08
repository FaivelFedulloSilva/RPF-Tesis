import pysam
from pysam import AlignmentFile
from dataclasses import dataclass


@dataclass
class BAMhandler:

    path: str
    bam_object: AlignmentFile = None

    def parse_bam(self, options):
        if options == 'b':
            self._bam_object = pysam.AlignmentFile(self.path, 'rb')
        elif options == 's':
            self._bam_object = pysam.AlignmentFile(self.path, 'rs')

    def region_coverage(self, contig: str, start: int, end: int):
        coverage = pysam.depth(self.path, '-a', '-r', f"{contig}:{start}-{end}")
        coverage = coverage.split('\n')
        coverage_return = []
        for line in coverage[:-1]:
            splitted_line = line.split('\t')
            coverage_return.append([splitted_line[1], splitted_line[2]])
        return coverage_return

# RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'
# TOTAL_PATH = './Data/totalRNA_01/accepted_hits_01.bam'

# contig = "chr15"
# start = 30903852
# end = 30903887
# regio = f"{contig}:{start}-{end}"
# bam = BAMhandler(RPF_PATH, 'b')
# print(bam.region_coverage(contig, start, end))
# counted = bam._bam_object.count(region=regio)
# print(counted)
# count_coverage = bam._bam_object.pileup(contig, start, end)
# for i in count_coverage:
#     print('->', i)

    

