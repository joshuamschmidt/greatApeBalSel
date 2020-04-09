# ncd_scores
paste the links below into the UCSC custom track config: https://genome.ucsc.edu/cgi-bin/hgCustom


track type=bigWig name="bonobo" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bonobo_ncd.log10p_ucsc.track.bw

track type=bigWig name="nigeria" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/ellioti_ncd.log10p_ucsc.track.bw

track type=bigWig name="western" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/verus_ncd.log10p_ucsc.track.bw

track type=bigWig name="eastern" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/schweinfurthii_ncd.log10p_ucsc.track.bw

track type=bigWig name="central" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/troglodytes_ncd.log10p_ucsc.track.bw

track type=bigWig name="gorilla" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/gorilla_ncd.log10p_ucsc.track.bw

track type=bigWig name="graueri" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/graueri_ncd.log10p_ucsc.track.bw

track type=bigWig name="abelii" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/abelii_ncd.log10p_ucsc.track.bw

track type=bigWig name="pygmaeus" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/pygmaeus_ncd.log10p_ucsc.track.bw


track type=bigWig name="bonobo" description="NCD" bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bonobo_ncd.log10p_ucsc.track.bw



coordinates are the central 100bp for each window.
score is -log10(pvalue)
** Note: for some reason converting from bedGraph to bigwig lost the default settings I made for the visualisation. 
User still has to manually change each track: full, 30 px, user defined range, but I like 0 - 2 ( p < 0.01 ) or 0-10. Otherwise the NCD values are indicated by the intenisty of shading, which is more difficult to parse visually.
