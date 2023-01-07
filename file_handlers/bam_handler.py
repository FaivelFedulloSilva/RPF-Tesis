import pysam


class BAMhandler:
    def __init__(self, path, options) -> None:
        self._path = path
    
        def _read_file(path, options):
            if options == 'b':
                return pysam.AlignmentFile(path, 'rb')
            elif options == 's':
                return pysam.AlignmentFile(path, 'rs')

        
        # self._bam_object = _read_file(path, options)

    def region_coverage(self, contig: str, start: int, end: int):
        coverage = pysam.depth(self._path, '-a', '-r', f"{contig}:{start}-{end}")
        coverage = coverage.split('\n')
        coverage_return = []
        for line in coverage[:-1]:
            splitted_line = line.split('\t')
            coverage_return.append([splitted_line[1], splitted_line[2]])
        return coverage_return

RPF_PATH = './Data/RPF_1/accepted_hits_01.bam'
TOTAL_PATH = './Data/totalRNA_01/accepted_hits_01.bam'


contig = "chr15"
start = 30903852
end = 30903887
regio = f"{contig}:{start}-{end}"
bam = BAMhandler(RPF_PATH, 'b')
print(bam.region_coverage(contig, start, end))
# counted = bam._bam_object.count(region=regio)
# print(counted)
# count_coverage = bam._bam_object.pileup(contig, start, end)
# for i in count_coverage:
#     print('->', i)

    

