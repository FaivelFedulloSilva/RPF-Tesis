import pysam


class BAMhandler:
    def __init__(self, path, options) -> None:
        self.__path = path
        
        def _read_file(path, options):
            if options == 'b':
                pysam.AlignmentFile(path, 'rb')
            elif options == 's':
                pysam.AlignmentFile(path, 'rs')

        
        self.__bam_object = _read_file(path, options)
    

    

