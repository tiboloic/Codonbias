# LT 23/06/2021
#
# extract synonymous variants SFS

import hail as hl
hl.init()
ht = hl.read_table('gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht')

# for call rate cutoff
AN_ADJ_FILTER = 0.8
def get_an_filter(ht):
	samples = {'female': 57787, 'male':67961}
	return (hl.case()
	.when(ht.locus.in_autosome_or_par(), ht.freq[0].AN >= AN_ADJ_FILTER * 2 * sum(samples.values()))
	.when(ht.locus.in_x_nonpar(), ht.freq[0].AN >= AN_ADJ_FILTER * (samples['male'] + samples['female'] * 2))
	.when(ht.locus.in_y_nonpar(), ht.freq[0].AN >= AN_ADJ_FILTER * samples['male'])
	.or_missing())

# keep only polymorphic synonymous variants in regions where the call rate is >80%
ht = ht.filter((ht.freq[0].AC > 0) & (ht.freq[0].AF < 1) & (ht.vep.most_severe_consequence=='synonymous_variant') & get_an_filter(ht))

# keep protein_coding transcripts where the variant is synonymous
ht = ht.annotate(tcs = ht.vep.transcript_consequences.filter(lambda x: (x.biotype=='protein_coding') & x.consequence_terms.contains('synonymous_variant')))

# get the set of codons
ht = ht.annotate(codons = hl.set(ht.tcs.map(lambda x: x.codons)))

# remove variants that have different codons in different transcripts (e.g. overlapping genes with different ORFs)
ht = ht.filter(ht.codons.size()==1)
ht = ht.select(codon = ht.codons.find(lambda x:True), AC=ht.freq[0].AC, AN=ht.freq[0].AN)

# bin allele counts per octave
ht = ht.annotate(cat = hl.floor(hl.log(ht.AC,2)))

# group by codon and summarize 
res = ht.group_by('codon','cat').aggregate(obs=hl.agg.count(), ansum=hl.agg.sum(ht.AN), anmax=hl.agg.max(ht.AN))

# save result
res.export('extract.synonymous.SFS.bycodon.tsv')
