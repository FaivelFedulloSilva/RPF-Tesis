import polars as pl
from polars import DataFrame

class GTFobject:

    def __init__(self, gtf: DataFrame) -> None:
        self._gtf = gtf


    def filter_by_feature(self, feature: str):
        """Return a new GTFobject that contains only the rows of type feature"""

        filtered_gtf = self._gtf.filter(
            pl.col("feature") == feature
        ) 
        return GTFobject(filtered_gtf)


    def get_transcripts_ids(self) -> list[str]:
        """Return a list of all the transcripts_id"""
        
        query = (self._gtf
                .lazy()
                .select(pl.col('transcript_id')).unique()
            )
        
        return query.collect().to_dict(as_series=False)['transcript_id']


    def get_transcripts_data(self, transcripts_ids: list[str]) -> DataFrame:
        query = (self._gtf
                .lazy()
                .groupby(by='transcript_id')
                .agg(
                    [pl.col('seqname').unique().first(), pl.col('strand'), pl.col('start'),  pl.col('end')]
                    )
                .filter(pl.col('transcript_id').is_in(transcripts_ids))
            ).collect()

        return query
    
    def get_transcript_id_data(self, transcript_id: str) -> DataFrame:
        return self.get_transcripts_data([transcript_id])

