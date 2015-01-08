
###################################################
### code chunk number 1: setup
###################################################
options(width=90)


###################################################
### code chunk number 2: biocLite (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("ShortRead", "VariantAnnotation")) # new packages
## biocLite()                                    # update packages


###################################################
### code chunk number 3: help-start (eval = FALSE)
###################################################
## help.start()


###################################################
### code chunk number 4: help (eval = FALSE)
###################################################
## library(ShortRead)
## ?readFastq


###################################################
### code chunk number 5: S4
###################################################
library(Biostrings)
showMethods(complement)


###################################################
### code chunk number 6: S4-showMethods (eval = FALSE)
###################################################
## showMethods(class="DNAStringSet", where=getNamespace("Biostrings"))


###################################################
### code chunk number 7: S4-help (eval = FALSE)
###################################################
## class ? DNAStringSet
## method ? "complement,DNAStringSet"


###################################################
### code chunk number 8: vignette (eval = FALSE)
###################################################
## vignette(package="useR2013")


###################################################
### code chunk number 9: scavenge (eval = FALSE)
###################################################
## ??readFastq
## library(Biostrings)
## ?alphabetFrequency
## class?GappedAlignments
## vignette(package="GenomicRanges")


### R code from vignette source 'SequencesAndRanges.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=90)
library(useR2013)


###################################################
### code chunk number 2: GRanges-genes
###################################################
genes = GRanges(seqnames=c("3R", "X"),
                 ranges=IRanges(
                     start=c(19967117, 18962306),
                     end=c(19973212, 18962925)),
                 strand=c("+", "-"),
                 seqlengths=c(`3R`=27905053L, `X`=22422827L))


###################################################
### code chunk number 3: GRanges-display
###################################################
genes


###################################################
### code chunk number 4: GRanges-help (eval = FALSE)
###################################################
## ?GRanges


###################################################
### code chunk number 5: GRanges-vignettes (eval = FALSE)
###################################################
## vignette(package="GenomicRanges")


###################################################
### code chunk number 6: ranges-ir
###################################################
ir &lt;- IRanges(start=c(7, 9, 12, 14, 22:24),
              end=c(15, 11, 12, 18, 26, 27, 28))


###################################################
### code chunk number 7: ranges-ir-plot
###################################################
png("ranges-ir-plot.png", width=800, height=160)
plotRanges(ir, xlim=c(5, 35), main="Original")
dev.off()
png("ranges-shift-ir-plot.png", width=800, height=160)
plotRanges(shift(ir, 5), xlim=c(5, 35), main="Shift")
dev.off()
png("ranges-reduce-ir-plot.png", width=800, height=160)
plotRanges(reduce(ir), xlim=c(5, 35), main="Reduce")
dev.off()


###################################################
### code chunk number 8: GRanges-mcols
###################################################
mcols(genes) &lt;- DataFrame(EntrezId=c("42865", "2768869"),
                          Symbol=c("kal-1", "CG34330"))


###################################################
### code chunk number 9: GRanges-metadata
###################################################
metadata(genes) &lt;-
    list(CreatedBy="A. User", Date=date())


###################################################
### code chunk number 10: GRangesList-eg-setup
###################################################
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb &lt;- TxDb.Dmelanogaster.UCSC.dm3.ensGene # alias
fbgn &lt;- exonsBy(txdb, "gene")["FBgn0039155"]
seqlevels(fbgn) &lt;- "chr3R"


###################################################
### code chunk number 11: GRangesList-eg
###################################################
fbgn


###################################################
### code chunk number 12: txdb
###################################################
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb &lt;- TxDb.Dmelanogaster.UCSC.dm3.ensGene # alias
ex0 &lt;- exonsBy(txdb, "gene")
head(table(elementLengths(ex0)))
ids &lt;- c("FBgn0002183", "FBgn0003360", "FBgn0025111", "FBgn0036449")
ex &lt;- reduce(ex0[ids])


###################################################
### code chunk number 13: gc-genome
###################################################
library(BSgenome.Dmelanogaster.UCSC.dm3)
nm &lt;- as.character(unique(seqnames(ex[[1]])))
chr &lt;- Dmelanogaster[[nm]]
v &lt;- Views(chr, start=start(ex[[1]]), end=end(ex[[1]]))


###################################################
### code chunk number 14: gcFunction-genome
###################################################
gcFunction(v)


###################################################
### code chunk number 15: gcFunction-definition
###################################################
gcFunction


### R code from vignette source 'ReadsAndAlignments.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=90)
library(useR2013)


###################################################
### code chunk number 2: fastq-format
###################################################
fls &lt;- dir(file.path(bigdata(), "fastq"), full=TRUE)
cat(noquote(readLines(fls[[1]], 4)), sep="\n")


###################################################
### code chunk number 3: ascii
###################################################
cat(rawToChar(as.raw(32+1:47)),
    rawToChar(as.raw(32+48:94)), sep="\n")


###################################################
### code chunk number 4: readFastq
###################################################
fastqDir &lt;- file.path(bigdata(), "fastq")
fastqFiles &lt;- dir(fastqDir, full=TRUE)
fq &lt;- readFastq(fastqFiles[1])
fq


###################################################
### code chunk number 5: sread
###################################################
head(sread(fq), 3)
head(quality(fq), 3)


###################################################
### code chunk number 6: width-ShortReadQ
###################################################
abc &lt;- alphabetByCycle(sread(fq))
abc[1:4, 1:8]


###################################################
### code chunk number 7: FastqSampler
###################################################
sampler &lt;- FastqSampler(fastqFiles[1], 1000000)
yield(sampler) # sample of 1000000 reads


###################################################
### code chunk number 8: qa (eval = FALSE)
###################################################
## ## Bioc 2.13 only; see ?qa for Bioc 2.12
## qas &lt;- qa(fastqFiles, type="fastq")
## rpt &lt;- report(qas, dest=tempfile())
## browseURL(rpt)


###################################################
### code chunk number 9: report (eval = FALSE)
###################################################
## rpt &lt;- system.file("GSM461176_81_qa_report", "index.html",
##                    package="useR2013")
## browseURL(rpt)


###################################################
### code chunk number 10: fastq-discovery
###################################################
dir(bigdata())
fls &lt;- dir(file.path(bigdata(), "fastq"), full=TRUE)


###################################################
### code chunk number 11: fastq-input-gc
###################################################
rm(fq); invisible(gc())


###################################################
### code chunk number 12: fastq-input
###################################################
fq &lt;- readFastq(fls[1])


###################################################
### code chunk number 13: gcC
###################################################
alf0 &lt;- alphabetFrequency(sread(fq), as.prob=TRUE, collapse=TRUE)
sum(alf0[c("G", "C")])


###################################################
### code chunk number 14: gc-reads
###################################################
gc &lt;- gcFunction(sread(fq))
hist(gc)


###################################################
### code chunk number 15: abc
###################################################
abc &lt;- alphabetByCycle(sread(fq))
matplot(t(abc[c("A", "C", "G", "T"),]), type="l")


###################################################
### code chunk number 16: abc-mclapply (eval = FALSE)
###################################################
## library(parallel)
## gc0 &lt;- mclapply(fls, function(fl) {
##   fq &lt;- readFastq(fl)
##   gc &lt;- gcFunction(sread(fq))
##   table(cut(gc, seq(0, 1, .05)))
## })
## ## simplify list of length 2 to 2-D array
## gc &lt;- simplify2array(gc0)
## matplot(gc, type="s")


###################################################
### code chunk number 17: SAM
###################################################
fl &lt;- system.file("extdata", "ex1.sam", package="Rsamtools")
strsplit(readLines(fl, 1), "\t")[[1]]


###################################################
### code chunk number 18: readGappedAlignments
###################################################
alnFile &lt;- system.file("extdata", "ex1.bam", package="Rsamtools")
aln &lt;- readGappedAlignments(alnFile)
head(aln, 3)


###################################################
### code chunk number 19: GappedAlignments-accessors
###################################################
table(strand(aln))
table(width(aln))
head(sort(table(cigar(aln)), decreasing=TRUE))


###################################################
### code chunk number 20: bam-ex-fls
###################################################
fls &lt;- dir(file.path(bigdata(), "bam"), ".bam$", full=TRUE)
names(fls) &lt;- sub("_.*", "", basename(fls))


###################################################
### code chunk number 21: bam-ex-input
###################################################
## input
aln &lt;- readGappedAlignments(fls[1])
xtabs(~seqnames + strand, as.data.frame(aln))


###################################################
### code chunk number 22: bam-ex-roi
###################################################
data(ex)             # from an earlier exercise


###################################################
### code chunk number 23: bam-ex-strand
###################################################
strand(aln) &lt;- "*"   # protocol not strand-aware


###################################################
### code chunk number 24: bam-ex-hits
###################################################
hits &lt;- findOverlaps(aln, ex)


###################################################
### code chunk number 25: qhits
###################################################
qhits &lt;- countQueryHits(hits)
table(qhits)


###################################################
### code chunk number 26: qhits-keep
###################################################
keep &lt;- which(qhits == 1)


###################################################
### code chunk number 27: bam-ex-cnt
###################################################
cnt &lt;- countSubjectHits(hits[queryHits(hits) %in% keep])


###################################################
### code chunk number 28: bam-count-fun
###################################################
counter &lt;-
    function(filePath, range)
{
    hits &lt;- findOverlaps(aln, ex)
    keep &lt;- which(countQueryHits(hits) == 1)
    cnts &lt;- countSubjectHits(hits[queryHits(hits) %in% keep])
    setNames(cnts, names(ex))
}


###################################################
### code chunk number 29: bam-count-all
###################################################
counts &lt;- sapply(fls, counter, ex)
counts


###################################################
### code chunk number 30: bam-count-mclapply (eval = FALSE)
###################################################
## if (require(parallel))
##     simplify2array(mclapply(fls, counter, ex))


###################################################
### code chunk number 31: gc-read
###################################################
param &lt;- ScanBamParam(what="seq")
seqs &lt;- scanBam(fls[[1]], param=param)
readGC &lt;- gcFunction(seqs[[1]][["seq"]])
hist(readGC)


### R code from vignette source 'RNASeq.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=90)
library(useR2013)


###################################################
### code chunk number 2: counts
###################################################
data(counts)
dim(counts)
grps &lt;- factor(sub("[1-4].*", "", colnames(counts)),
               levels=c("untreated", "treated"))
pairs &lt;- factor(c("single", "paired", "paired",
                  "single", "single", "paired", "paired"))
pData &lt;- data.frame(Group=grps, PairType=pairs,
                    row.names=colnames(counts))


###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)
dge &lt;- DGEList(counts, group=pData$Group)
dge &lt;- calcNormFactors(dge)


###################################################
### code chunk number 4: DEGList-filter
###################################################
m &lt;- sweep(dge$counts, 2, 1e6 / dge$samples$lib.size, `*`)
ridx &lt;- rowSums(m &gt; 1) &gt;= 2
table(ridx)                        # number filtered / retained
dge &lt;- dge[ridx,]


###################################################
### code chunk number 5: design
###################################################
(design &lt;- model.matrix(~ Group, pData))


###################################################
### code chunk number 6: common.dispersion
###################################################
dge &lt;- estimateTagwiseDisp(dge)
mean(sqrt(dge$tagwise.dispersion))


###################################################
### code chunk number 7: glmFit
###################################################
fit &lt;- glmFit(dge, design)


###################################################
### code chunk number 8: lrt
###################################################
lrTest &lt;- glmLRT(fit, coef=2)


###################################################
### code chunk number 9: topTags
###################################################
tt &lt;- topTags(lrTest, n=10)
tt[1:3,]


###################################################
### code chunk number 10: sanity
###################################################
sapply(rownames(tt$table)[1:4],
       function(x) tapply(counts[x,], pData$Group, mean))


### R code from vignette source 'ChIPSeq.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=90)
library(useR2013)


###################################################
### code chunk number 2: chipseq-report
###################################################
rpt &lt;- system.file("GSE30263_qa_report", "index.html",
           package="useR2013", mustWork=TRUE)
if (interactive())
    browseURL(rpt)


###################################################
### code chunk number 3: chipseq-halfpeak-stamFile
###################################################
stamFile &lt;- system.file("data", "stam.Rda", package="useR2013")
load(stamFile)
stam


###################################################
### code chunk number 4: chipseq-stam
###################################################
head(colData(stam), 3)
head(rowData(stam), 3)
xtabs(~Replicate + CellLine, colData(stam))[,1:5]


###################################################
### code chunk number 5: chipseq-stam-detect
###################################################
m &lt;- assays(stam)[["Tags"]] &gt; 0 # peaks detected...
peaksPerSample &lt;- table(rowSums(m))
head(peaksPerSample)
tail(peaksPerSample)


###################################################
### code chunk number 6: chipseq-stam-similarity-1
###################################################
library(bioDist)                   # for cor.dist
m &lt;- asinh(assays(stam)[["Tags"]]) # transformed tag counts
d &lt;- cor.dist(t(m))                # correlation distance
h &lt;- hclust(d)                     # hierarchical clustering


###################################################
### code chunk number 7: chipseq-stam-similarity-plot (eval = FALSE)
###################################################
## plot(h, cex=.8, ann=FALSE)


###################################################
### code chunk number 8: chipseq-stam-similarity
###################################################
png("chipseq-stam-similarity.png", width=1280)
opar &lt;- par(mar=c(0, 0, 0, 0))
plot(h, axes=FALSE, ann=FALSE)
par(opar)
invisible(dev.off())


###################################################
### code chunk number 9: chipseq-CTCF-PWM-setup
###################################################
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqLogo)
library(lattice)
pwm &lt;- getJASPAR("MA0139.1") # useR2013::getJASPAR


###################################################
### code chunk number 10: chipseq-CTCF-PWM-binding
###################################################
chrid &lt;- "chr1"
hits &lt;-matchPWM(pwm, Hsapiens[[chrid]]) # '+' strand
scores &lt;- PWMscoreStartingAt(pwm, subject(hits), start(hits))


###################################################
### code chunk number 11: chipseq-CTCF-PWM-densityplot (eval = FALSE)
###################################################
## densityplot(scores, xlim=range(scores), pch="|")


###################################################
### code chunk number 12: CTCF-PWM-found
###################################################
cm &lt;- consensusMatrix(hits)[1:4,]
seqLogo(makePWM(scale(cm, FALSE, colSums(cm))))


### R code from vignette source 'Annotation.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=90)
library(useR2013)


###################################################
### code chunk number 2: select
###################################################
cols(org.Dm.eg.db)
keytypes(org.Dm.eg.db)
uniprotKeys &lt;- head(keys(org.Dm.eg.db, keytype="UNIPROT"))
cols &lt;- c("SYMBOL", "PATH")
select(org.Dm.eg.db, keys=uniprotKeys, cols=cols, keytype="UNIPROT")


###################################################
### code chunk number 3: select-kegg
###################################################
kegg &lt;- select(org.Dm.eg.db, "00310", c("UNIPROT", "SYMBOL"), "PATH")
nrow(kegg)
head(kegg, 3)


###################################################
### code chunk number 4: chipseq-anno-data
###################################################
stamFile &lt;- system.file("data", "stam.Rda", package="useR2013")
load(stamFile)


###################################################
### code chunk number 5: chipseq-anno-common
###################################################
ridx &lt;- rowSums(assays(stam)[["Tags"]] &gt; 0) == ncol(stam)
peak &lt;- rowData(stam)[ridx]


###################################################
### code chunk number 6: chipseq-anno-centers
###################################################
peak &lt;- resize(peak, width=1, fix="center")


###################################################
### code chunk number 7: chipseq-anno-tss
###################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
tx &lt;- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
tss &lt;- resize(tx, width=1)


###################################################
### code chunk number 8: chipseq-anno-tss-dist
###################################################
idx &lt;- nearest(peak, tss)
sgn &lt;- as.integer(ifelse(strand(tss)[idx] == "+", 1, -1))
dist &lt;- (start(peak) - start(tss)[idx]) * sgn


###################################################
### code chunk number 9: chipseq-anno-tss-dist
###################################################
bound &lt;- 1000
ok &lt;- abs(dist) &lt; bound
dist &lt;- dist[ok]
table(sign(dist))


###################################################
### code chunk number 10: anno-tss-disthist
###################################################
griddensityplot &lt;-
    function(...)
    ## 'panel' function to plot a grid underneath density
{
    panel.grid()
    panel.densityplot(...)
}
print(densityplot(dist[ok], plot.points=FALSE,
    panel=griddensityplot,
    xlab="Distance to Nearest TSS"))


###################################################
### code chunk number 11: tss-dist-func
###################################################
distToTss &lt;-
    function(peak, tx)
{
    peak &lt;- resize(peak, width=1, fix="center")
    tss &lt;- resize(tx, width=1)
    idx &lt;- nearest(peak, tss)
    sgn &lt;- as.numeric(ifelse(strand(tss)[idx] == "+", 1, -1))
    (start(peak) - start(tss)[idx]) * sgn
}


###################################################
### code chunk number 12: chipseq-anno-seq
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)
ridx &lt;- rowSums(assays(stam)[["Tags"]] &gt; 0) == ncol(stam)
ridx &lt;- ridx &amp; (seqnames(rowData(stam)) == "chr6")
pk6 &lt;- rowData(stam)[ridx]
seqs &lt;- getSeq(Hsapiens, pk6)
head(seqs, 3)


###################################################
### code chunk number 13: chipseq-tss-dist-2
###################################################
pwm &lt;- useR2013::getJASPAR("MA0139.1")
hits &lt;- lapply(seqs, matchPWM, pwm=pwm)
hasPwmMatch &lt;- sapply(hits, length) &gt; 0
dist &lt;- distToTss(pk6, tx)

ok &lt;- abs(dist) &lt; bound
df &lt;- data.frame(Distance = dist[ok], HasPwmMatch = hasPwmMatch[ok])
print(densityplot(~Distance, group=HasPwmMatch, df,
    plot.points=FALSE, panel=griddensityplot,
    auto.key=list(
      columns=2,
      title="Has Position Weight Matrix?",
      cex.title=1),
    xlab="Distance to Nearest Tss"))


###################################################
### code chunk number 14: readVcf
###################################################
library(VariantAnnotation)
fl &lt;- system.file("extdata", "chr22.vcf.gz",
                  package="VariantAnnotation")
(hdr &lt;- scanVcfHeader(fl))
info(hdr)[c("VT", "RSQ"),]


###################################################
### code chunk number 15: readVcf
###################################################
(vcf &lt;- readVcf(fl, "hg19"))
head(rowData(vcf), 3)


###################################################
### code chunk number 16: renameSeqlevels
###################################################
rowData(vcf) &lt;- renameSeqlevels(rowData(vcf), c("22"="ch22"))


###################################################
### code chunk number 17: dbSNP
###################################################
library(SNPlocs.Hsapiens.dbSNP.20101109)
snpFilt &lt;- useR2013::dbSNPFilter("SNPlocs.Hsapiens.dbSNP.20101109")
inDbSNP &lt;- snpFilt(vcf)
table(inDbSNP)


###################################################
### code chunk number 18: SNP-quality
###################################################
metrics &lt;-
    data.frame(inDbSNP=inDbSNP, RSQ=info(vcf)$RSQ)


###################################################
### code chunk number 19: RSQ-plot
###################################################
library(ggplot2)
ggplot(metrics, aes(RSQ, fill=inDbSNP)) +
    geom_density(alpha=0.5) +
    scale_x_continuous(name="MaCH / Thunder Imputation Quality") +
    scale_y_continuous(name="Density") +
    theme(legend.position="top")


</pre><embed id="xunlei_com_thunder_helper_plugin_d462f475-c18e-46be-bd10-327458d045bd" type="application/thunder_download_plugin" height="0" width="0"></body></html>Ztext/plain    ( F ] l ~ ? ? ? ?a?             
              a?