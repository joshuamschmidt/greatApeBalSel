# ncd_scores
paste the links below into the UCSC custom track config: https://genome.ucsc.edu/cgi-bin/hgCustom


track type=bigWig name="bonobo" description="NCD_bonobo" color="0,0,128" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/bonobo_ncd.log10p_ucsc.track.bw

track type=bigWig name="nigeria" description="NCD_NC_chimp" color="123,104,238" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/ellioti_ncd.log10p_ucsc.track.bw

track type=bigWig name="eastern" description="NCD_eastern_chimp" color="70,130,180" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/schweinfurthii_ncd.log10p_ucsc.track.bw

track type=bigWig name="central" description="NCD_central_chimp" color="30,144,255" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/troglodytes_ncd.log10p_ucsc.track.bw

track type=bigWig name="western" description="NCD_western_chimp" color="72,61,139" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/verus_ncd.log10p_ucsc.track.bw

track type=bigWig name="gorilla" description="NCD_western_llg" color="255,69,0" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/gorilla_ncd.log10p_ucsc.track.bw

track type=bigWig name="graueri" description="NCD_eastern_llg" color="255,127,80" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/graueri_ncd.log10p_ucsc.track.bw

track type=bigWig name="abelii" description="NCD_sumatran_orang" color="139,0,0" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/abelii_ncd.log10p_ucsc.track.bw

track type=bigWig name="pygmaeus" description="NCD_bornean_orang" color="220,20,60" autoScale=off alwaysZero=on visibility=2 viewLimits=0:10 yLineMark=3.3 yLineOnOff=on maxHeightPixels=30:30:30 bigDataUrl=https://raw.githubusercontent.com/joshuamschmidt/ncd_scores/master/bigWig/pygmaeus_ncd.log10p_ucsc.track.bw


coordinates are the central 100bp for each window.
score is -log10(pvalue)
