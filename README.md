# ncd_scores
paste this into the UCSC custom track config: https://genome.ucsc.edu/cgi-bin/hgCustom

https://github.com/joshuamschmidt/ncd_scores/raw/master/bonobo_ncd.log10p_ucsc.track.bw

https://github.com/joshuamschmidt/ncd_scores/raw/master/ellioti_ncd.log10p_ucsc.track.bw

https://github.com/joshuamschmidt/ncd_scores/raw/master/verus_ncd.log10p_ucsc.track.bw

https://github.com/joshuamschmidt/ncd_scores/raw/master/schweinfurthii_ncd.log10p_ucsc.track.bw

https://github.com/joshuamschmidt/ncd_scores/raw/master/troglodytes_ncd.log10p_ucsc.track.bw

https://github.com/joshuamschmidt/ncd_scores/raw/master/gorilla_ncd.log10p_ucsc.track.bw

https://github.com/joshuamschmidt/ncd_scores/raw/master/graueri_ncd.log10p_ucsc.track.bw

https://github.com/joshuamschmidt/ncd_scores/raw/master/abelii_ncd.log10p_ucsc.track.bw

https://github.com/joshuamschmidt/ncd_scores/raw/master/pygmaeus_ncd.log10p_ucsc.track.bw

coordinates are the central 100bp for each window.
score is -log10(pvalue)
** Note: for some reason converting from bedGraph to bigwig lost the default settings I made for the visualisation. 
User still has to manually change each track: full, 30 px, user defined range, but I like 0 - 2 ( p < 0.01 ) or 0-10.
