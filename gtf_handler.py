from gtfparse import read_gtf 
import polars as pl

class GTFhandler:
    def __init__(self, path: str) -> None:
        self.__path: str = path
        self.__gtf_object = read_gtf(path)

    def filter_by_feature(self, feature: str):
        filtered_gtf = self.__gtf_object.filter(
            pl.col("feature") == feature
        )
        return filtered_gtf

gtf = GTFhandler('./Data/genesFiltrada.gtf')
# print(gtf._gtf_object.describe())
gtf.filter_by_feature('CDS')     
