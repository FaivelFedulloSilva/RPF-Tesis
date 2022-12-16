from gtfparse import read_gtf 

class GTFhandler:
    def __init__(self, path) -> None:
        self.__path = path
        self.__gtf_object = read_gtf(path)
        print(self.__gtf_object)

    def filter_by_feature(self, feature: str):
        print(self.__gtf_object[self.__gtf_object.feature == feature])
        print(self.__gtf_object.feature.unique())

gtf = GTFhandler('./Data/genesFiltrada.gtf')
gtf.filter_by_feature('CDS')     
