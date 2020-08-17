# ncd_scores
paste the links below into the UCSC custom track config: https://genome.ucsc.edu/cgi-bin/hgCustom


track type=bigWig name="bonobo" description="NCD_bonobo" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/bonobo_ncd.log10p_ucsc.track.bw

track type=bigWig name="nigeria" description="NCD_NC_chimp" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/ellioti_ncd.log10p_ucsc.track.bw

track type=bigWig name="western" description="NCD_western_chimp" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/verus_ncd.log10p_ucsc.track.bw

track type=bigWig name="eastern" description="NCD_eastern_chimp" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/schweinfurthii_ncd.log10p_ucsc.track.bw

track type=bigWig name="central" description="NCD_central_chimp" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/troglodytes_ncd.log10p_ucsc.track.bw

track type=bigWig name="gorilla" description="NCD_western_llg" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/gorilla_ncd.log10p_ucsc.track.bw

track type=bigWig name="graueri" description="NCD_eastern_llg" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/graueri_ncd.log10p_ucsc.track.bw

track type=bigWig name="abelii" description="NCD_sumatran_orang" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/abelii_ncd.log10p_ucsc.track.bw

track type=bigWig name="pygmaeus" description="NCD_bornean_orang" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/pygmaeus_ncd.log10p_ucsc.track.bw


coordinates are the central 100bp for each window.
score is -log10(pvalue)
** Note: for some reason converting from bedGraph to bigwig lost the default settings I made for the visualisation. 
User still has to manually change each track: full, 30 px, user defined range, but I like 0 - 2 ( p < 0.01 ) or 0-10. Otherwise the NCD values are indicated by the intenisty of shading, which is more difficult to parse visually.
