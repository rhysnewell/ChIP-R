"""
This python module reads in sam files from RNA-seq experiment and processes
them and RNA-seq data1
"""

from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import itertools
import operator
import math
from scipy import stats
from numpy import array, empty
import scipy.cluster.hierarchy as sch


def sam_reader(filename):
    """Mandatory fields are QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL
for more info http://samtools.github.io/hts-specs/SAMv1.pdf """
    data=[]
    f= open(filename,'r')
    for row in f:
        if row.startswith('@'): # skip the header
            pass
        else:
            info=row.strip().split('\t')
            data.append(info)
    return data


def base_percentages(reads):
    "reports base percentage  %A,%T,%C,%G "
    all_seqs=[]
    for read in reads:
        seq=read[9]
        seq=[seq[i:i+1] for i in range(0,len(seq),1)]
        for nuc in seq:
            all_seqs.append(nuc)
    counts=dict(Counter(all_seqs))
    nucs=list(counts.keys())
    freqs={}
    for nuc in nucs:
        freqs[nuc]=float(counts[nuc])/sum(counts.values())
    return freqs


def numberofreads(reads):
    """Incremented for every sequence-containing line in the sam file, regardless of whether it represents an alignment.
for some files, this is not actually the number of reads. indeed, this may be a poor name for this stat"""
    return len(reads)


def mapped_reads(reads,paired_end=True):
    """If duplicate tracking was enabled via -D, then this attempts to recapitulate the number of unique, mapped, probe-id's in the original sam file. It is multiplied by 2 for paired-end data with duplicate read id's.
The idea is that if you divide this by the number of reads in the fastq you aligned (possibly from the output of fastq-stats),
you will get an accurate "percentage of reads aligned" statistic.
"mapped" is something with a non-negative position, and a "non-asterisk" cigar string."""
    mapped_reads=[]
    store_reads=[]
    for read in reads:
        if read[3]>0 and read[5]!='*':
            mapped_reads.append(read[0])
            store_reads.append(read)
    mapped=set(mapped_reads)
    list_mapped=list(mapped)
    if paired_end==True:
        mapped=len(mapped)+len(mapped)
    else:
        mapped=len(mapped)
    print("number of mapped reads",mapped)
    return store_reads


def mappedBases(mapped_reads):
    """Total number of mapped bases in sam file"""
    seq=""
    for read in mapped_reads:
        seq=seq+read[9]
    return len(seq)

def forward(mapped_reads):
    """The number of lines in the sam file that were aligned to the "forward" strand. No accounting is done on duplicates."""
    forward=[read for read in mapped_reads if read[9]>0]
    return forward


def reverse(mapped_reads):
    """The number of lines in the sam file that were aligned to the "reverse" strand. No accounting is done on duplicates."""
    reverse=[read for read in mapped_reads if read[9]<0]
    return reverse



########Qualities and STATS


def subgroups(mapped_reads):
    """form groups p<1e-3 one group,1e-3<=p<1e-2 one group,1e-2<=p<1 one group a total of three groups"""
    group1=[]
    group2=[]
    group3=[]
    for read in mapped_reads:
        if int(read[4])>29:
            group1.append(read)
        elif int(read[4])<=29 and int(read[4])>17:
            group2.append(read)
        elif int(read[4])<=17:
            group3.append(read)
        else:
            pass
    print(len(group1),"in p<1e-3 group")
    print(len(group2),"in 1e-3<=p<1e-2 group")
    print(len(group3),"in 1e-2<=p<1 group")
    return group1,group2,group3



def dinuc_freq(mapped_reads):
    "reports dinucleotide composition using p(Rho) statistics for overrepresentation"
    all_seqs=[]
    for read in mapped_reads:
        seq=read[9]
        seq=[seq[i:i+1] for i in range(0,len(seq),1)]
        for nuc in seq:
            all_seqs.append(nuc)
    counts=dict(Counter(all_seqs))
    nucs=list(counts.keys())
    freqs={}
    for nuc in nucs:
        freqs[nuc]=float(counts[nuc])/sum(counts.values())
    all_seqs=[]
    for read in mapped_reads:
        seq=read[9]
        seq=[seq[i:i+2] for i in range(0,len(seq),2)]
        for nuc in seq:
            all_seqs.append(nuc)
    counts=dict(Counter(all_seqs))
    dinucs=list(counts.keys())
    dinuc_counts={}
    for i in dinucs:
        val=float(counts[i])/sum(counts.values())
        dinuc_counts[i]=val/(freqs[i[0]]*freqs[i[1]]) # p-values
    return dinuc_counts


def PercentReadsAligned(group1,group2,group3,numfastq):
    """Provide a list of mapped_reads and the number of reads in the fastq file"""
    mapped_reads=group1+group2+group3
    Mapped=len(mapped_reads)/float(numfastq)
    Unmapped=1-float(Mapped)
##    print "Mapping stats"
##    print"p<1e-3", len(group1)/float(numfastq)
##    print"1e-3<=p<1e-2",len(group2)/float(numfastq)
##    print "1e-2<=p<1",len(group3)/float(numfastq)
##    print "Unmapped",Unmapped
    labels="p<1e-3","1e-3<=p<1e-2","1e-2<=p<1","Unmapped"
    x=[len(group1)/float(numfastq),len(group2)/float(numfastq),len(group3)/float(numfastq),Unmapped]
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.pie(x,labels=labels,autopct='%1.1f%%', shadow=True)
    plt.title('Mapping stats')
    plt.show()
    return Mapped



def length_stats(group1,group2,group3):
    """returns basic stats relating to the lengths of the reads
Calculations are based on the the length of the (possibly hard-clipped) sequence in the sam file."""
    reads=[group1,group2,group3]
    data=[]
    for i in range(0,len(reads)):
        lengths=[]
        for read in reads[i]:
            if int(read[8])<0:
                length=-1*int(read[8])
            else:
                length=int(read[8])
            lengths.append(length)
        mean_len=np.mean(lengths)
        print("group"+str(i+1)+"mean",mean_len)
        max_len=np.max(lengths)
        print("group"+str(i+1)+"max length",max_len)
        min_len=np.min(lengths)
        print("group"+str(i+1)+"min length",min_len)
        data.append(["group"+str(i+1),mean_len,max_len,min_len])
    return data

def plot_length_distrib(group,name):
    """distribution of lengths of all the sam reads"""
    lengths=[]
    for read in group:
        if int(read[8])<0:
            length=-1*int(read[8])
        else:
            length=int(read[8])
        lengths.append(length)
    ##Visualize length distribution
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    n, bins, patches = plt.hist(lengths,100, normed=0, facecolor='g')
    plt.xlabel("lengths")
    plt.ylabel("number of mapped reads")
    plt.title(name)
    plt.show()

def inv_logit(p):
    return 10**(p/-10)



def plot_base_composition(reads,sym):
    "reports nucelotide frequencies at each position in the sam sequences"
    #DNA_Alphabet=["A","C","T","G","N"]
    all_nucs=[]
    for read in reads:
        nucs={}#dictionary to store nucleotide data
        seq=read[9]
        for i in range(0,len(seq)):
            nucs[str(i+1)]=seq[i]
        all_nucs.append(nucs)
    all_items=[]
    counts=[]
    pos=list(range(1,len(seq)+1))
    for dicts in all_nucs:
        for item in list(dicts.items()):
            all_items.append(item)
    all_items.sort(key=operator.itemgetter(0))
    groups= [list(map(operator.itemgetter(1),list(group))) for key, group in itertools.groupby(all_items, operator.itemgetter(0))]
    for group in groups:
        counts.append(group.count(sym))
    print(counts)
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.bar(pos,counts,facecolor='g')
    plt.xlabel("Position")
    plt.ylabel("number of mapped reads")
    plt.title(sym)
    plt.show()
    return counts

#####################################################
#Transcript reader

def raw_count_reader(filename):
    """Count the raw counts in the file"""
    data={}
    f= open(filename,'r')
    for row in f:
        if row.startswith('t1'): # skip the header
            pass
        else:
            info=row.strip().split('\t')
            data[info[0]]=[int(info[1]),int(info[2]),int(info[3]),int(info[4]),float(info[5])] #t1,rept1,t10,rept10,len
    return data


#####Normalisation methods


def get_RPKM(data,num_map1,num_map2,num_map3,num_map4):
    """provide number of mapped reads for the two groups of interest and raw count data .This method provides length normalisation to prevent length and total count bias"""
    all_rpkms=[];final={}
    for i,s,ii,ss,v in list(data.values()):
        rpkms=[]
        num_mapped_reads=[num_map1,num_map2,num_map3,num_map4]
        vals=[i,s,ii,ss]
        lengths=[v,v,v,v]
        for n in range(0,len(vals)):
            if vals[n]==0:
                rpkm=0
                rpkms.append(rpkm)
            else:
                #perform RPKM calc
                rpkm= float(vals[n])/(lengths[n]*(float(num_mapped_reads[n])/10**6))
                rpkms.append(rpkm)
        all_rpkms.append(rpkms)
    #return gene names and rpkms
    for i in range(0,len(list(data.keys()))):
        final[list(data.keys())[i]]=[float(all_rpkms[i][0]),float(all_rpkms[i][1]),float(all_rpkms[i][2]),float(all_rpkms[i][3])]
    return final

def write_RPKM_data(RPKM_data,filename):
    """write RPKM data to a file"""
    f=open(filename,'w')
    for i in range(0,len(RPKM_data)):
        f.write("%s\t%d\t%d\t%d\t%d\n"%(list(RPKM_data.keys())[i],int(list(RPKM_data.values())[i][0]),int(list(RPKM_data.values())[i][1]),int(list(RPKM_data.values())[i][2]),int(list(RPKM_data.values())[i][3])))
    f.close()




###############Visualize replicates to determine degree of biological variation

def pearson_def(x, y):
    """Pearson correlation coefficient R"""
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = np.mean(x)
    avg_y = np.mean(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / math.sqrt(xdiff2 * ydiff2)


def plotreprpkm(rpkm_data,timepoint):
    """plot showing level of agreement between technical replicates for RPKM between replicates and plots coefficient of determination"""
    one=[]
    two=[]
    if timepoint=="t1":
        for i in range(0,len(list(rpkm_data.values()))):
            one.append(int(list(rpkm_data.values())[i][0]))
            two.append(int(list(rpkm_data.values())[i][1]))
    else:
        for i in range(0,len(list(rpkm_data.values()))):
            one.append(int(list(rpkm_data.values())[i][2]))
            two.append(int(list(rpkm_data.values())[i][3]))
    plt.plot(one,two,'o')
    pcc=pearson_def(one,two)
    R2=pcc**2
    name="""Technical Replicates
R2="""+str(R2)
    m,b= np.polyfit(one,two,1)
    plt.figure(1, figsize=(8,8))
    plt.plot(one, np.array(one)*m +b,'r-')
    plt.text(3000, max(two)-1000,name , fontsize=12)
    plt.xlabel("RPKM replicate 1")
    plt.ylabel("RPKM replicate 2")
    plt.title(timepoint)
    plt.show()


def plotMAreprpkm(rpkm_data,timepoint):
    """MA Plot of log(RPKM) vs Average log(RPKM) of replicates"""
    m=[]
    a=[]
    if timepoint=="t1":
        for i in range(0,len(list(rpkm_data.values()))):
            y=np.log2(list(rpkm_data.values())[i][0]+1)-np.log2(list(rpkm_data.values())[i][1]+1)
            x=(np.log2(list(rpkm_data.values())[i][0]+1)+np.log2(list(rpkm_data.values())[i][1]+1))/2
            m.append(y)
            a.append(x)
    else:
        for i in range(0,len(list(rpkm_data.values()))):
            y=np.log2(list(rpkm_data.values())[i][2]+1)-np.log2(list(rpkm_data.values())[i][3]+1)
            x=(np.log2(list(rpkm_data.values())[i][2]+1)+np.log2(list(rpkm_data.values())[i][3]+1))/2
            m.append(y)
            a.append(x)
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.plot(a,m,'o')
    plt.axhline(np.mean(m)+1.96*np.std(m),color="green",label="avg diff +1.96(std diff)")
    plt.axhline(np.mean(m)-1.96*np.std(m),color="green",label="avg diff -1.96(std diff)")
    plt.xlabel("Average log(RPKM) of replicates")
    plt.ylabel("Difference in log(RPKM) of replicates")
    plt.legend(loc="lower right")
    plt.title(timepoint)
    plt.show()



def get_cv(data1,condition):
    cvs=[]
    if condition=="t1":
        for i in range(0,len(list(data1.values()))):
            mean = np.mean([list(data1.values())[i][0],list(data1.values())[i][1]])
            std=np.std([list(data1.values())[i][0],list(data1.values())[i][1]])
            if mean==0.0 and std==0.0:
                pass
            else:
                cv=float(mean+1)/(std+1)
                cvs.append(cv)
    else:
        for i in range(0,len(list(data1.values()))):
            mean = np.mean([list(data1.values())[i][2],list(data1.values())[i][3]])
            std=np.std([list(data1.values())[i][2],list(data1.values())[i][3]])
            if mean==0.0 and std==0.0:
                pass
            else:
                cv=float(mean+1)/(std+1)
                cvs.append(cv)
    return cvs


def get_boxplots(norm,original):
    """distribution of the coeficient of variation across samples (replicates) normalised using the methods provided"""
    bp=plt.boxplot([norm,original],notch=False, patch_artist=True)
    for box in bp['boxes']:
        box.set(color="red")
        box.set(color="blue")
    plt.ylabel("coefficient of variation")
    plt.xlabel("Methods")
    my_xticks = ['RPKM','raw counts']
    x=[1,2]
    plt.xticks(x,my_xticks)
    plt.ylim(0,400)
    plt.show()


def plotavg_cv(norm,original):
    """distribution of the coeficient of variation across samples (replicates) normalised using the methods provided"""
    x=[1,2]
    y=[np.mean(norm),np.mean(original)]
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.bar(x[0],y[0],color="red",label="RPKM")
    plt.bar(x[1],y[1],color="blue",label="Raw counts")
    plt.ylabel("Average coefficient of variation")
    plt.xlabel("Methods")
    ax.xaxis.set_ticklabels([])
    plt.legend(loc="upper right")
    plt.show()


def plotMA(rpkm_data,cutoff=[-1.5,1.5]):
    """Produce MA plot using logfold as cutoff"""
    logfc=[]
    avg_rpkm=[]
    sig_logfc=[]
    sig_avg_rpkm=[]
    logfc2=[]
    avg_rpkm2=[]
    sig_logfc2=[]
    sig_avg_rpkm2=[]
    for i,ii,s,ss in list(rpkm_data.values()):
        fc=np.log2(float(s+1)/(i+1))
        if fc<cutoff[0] or fc>cutoff[1]:
            sig_logfc.append(fc)
            sig_avg_rpkm.append(np.log2(s+1)+np.log2(i+1)/2)
        else:
            logfc.append(fc)
            avg_rpkm.append(np.log2(s+1)+np.log2(i+1)/2)
    for i,ii,s,ss in list(rpkm_data.values()):
        fc2=np.log2(float(ss+1)/(ii+1))
        if fc2<cutoff[0] or fc2>cutoff[1]:
            sig_logfc2.append(fc2)
            sig_avg_rpkm2.append(np.log2(ss+1)+np.log2(ii+1)/2)
        else:
            logfc2.append(fc2)
            avg_rpkm2.append(np.log2(ss+1)+np.log2(ii+1)/2)
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.plot(avg_rpkm,logfc,'o',color="blue",label="rep1")
    plt.plot(avg_rpkm2,logfc2,'x',color="blue",label="rep2")
    plt.plot(sig_avg_rpkm,sig_logfc,'o',color="red",label="sig rep1")
    plt.plot(sig_avg_rpkm2,sig_logfc2,'x',color="red",label="sig rep2")
    plt.axhline(cutoff[0],color="orange")
    plt.axhline(cutoff[1],color="orange")
    plt.ylabel("Fold Change (log2)")
    plt.xlabel("Average RPKM (log2)")
    plt.title("MA plot")
    plt.legend(loc="upper left")
    plt.show()

def plotMA_pval(rpkm_data,cutoff=0.05):
    """Produce MA plot using the pvalue as cutoff"""
    logfc=[]
    avg_rpkm=[]
    sig_logfc=[]
    sig_avg_rpkm=[]
    logfc2=[]
    avg_rpkm2=[]
    sig_logfc2=[]
    sig_avg_rpkm2=[]
    for i,ii,s,ss,pval in list(rpkm_data.values()):
        fc=np.log2(float(s+1)/(i+1))
        if float(pval)<cutoff:
            sig_logfc.append(fc)
            sig_avg_rpkm.append(np.log2(s+1)+np.log2(i+1)/2)
        else:
            logfc.append(fc)
            avg_rpkm.append(np.log2(s+1)+np.log2(i+1)/2)
    for i,ii,s,ss,pval in list(rpkm_data.values()):
        fc2=np.log2(float(ss+1)/(ii+1))
        if float(pval)<cutoff:
            sig_logfc2.append(fc2)
            sig_avg_rpkm2.append(np.log2(ss+1)+np.log2(ii+1)/2)
        else:
            logfc2.append(fc2)
            avg_rpkm2.append(np.log2(ss+1)+np.log2(ii+1)/2)
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.plot(avg_rpkm,logfc,'o',color="blue",label="rep1")
    plt.plot(avg_rpkm2,logfc2,'o',color="blue",label="rep2")
    plt.plot(sig_avg_rpkm,sig_logfc,'o',color="red",label="sig rep1")
    plt.plot(sig_avg_rpkm2,sig_logfc2,'x',color="red",label="sig rep2")
    plt.ylabel("Fold Change (log2)")
    plt.xlabel("Average RPKM (log2)")
    plt.title("MA plot")
    plt.legend(loc="upper left")
    plt.show()


#####DE expression statistical test (T-Test, ANOVA and FDR)


def Welcht(rpkm):
    """Performs Welchs T-statistic (one-tailed)"""
    ts=[]
    result={}
    for i,ii,s,ss in list(rpkm.values()):
        sd1=np.std([i,ii])
        sd2=np.std([s,ss])
        t=(np.mean([s,ss])-np.mean([i,ii]))/(math.sqrt(((float(sd2)/2)+(float(sd1)/2))))
        ts.append(t)
    pvals=[]
    for t in ts:
        pval = stats.t.sf(np.abs(t), 2-1)
        if pval==float('nan'):
            pval=1
            pvals.append(pval)
        else:
            pval=pval
            pvals.append(pval)
    corr_pvals=correct_pvalues_for_multiple_testing(pvals, correction_type = "Benjamini-Hochberg")
    for i in range(0,len(list(rpkm.values()))):
        result[list(rpkm.keys())[i]]=[list(rpkm.values())[i][0],list(rpkm.values())[i][1],list(rpkm.values())[i][2],list(rpkm.values())[i][3],corr_pvals[i]]
    return result



def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues


####Method Run hiearachical clustering on the correlation matrix (of differentially expressed genes) -Coexpression

def cluster_data(data_matrix,genenames,timepoint):
    "One replicates at a specific time point"
    D = np.zeros([np.shape(data_matrix)[0],1])
    ##generate a distance matrix
    for i in range(np.shape(data_matrix)[0]):
        for j in range(1):
            D[i,j] = abs(data_matrix[i] - data_matrix[j])**2  #use Wards method (other methods could be implemented here)
    labels=list('' for i in range(np.shape(data_matrix)[0]))
    for i in range(np.shape(data_matrix)[0]):
        labels[i]=str(i)+","+str(genenames[i])
    fig=plt.figure(1, figsize=(17,8))
    linked = sch.linkage(D, method='centroid')
    dend = sch.dendrogram(linked, orientation='right',labels=labels) # sets the oirentation root at the right
    plt.title(timepoint)
    fig.savefig(timepoint+'dendogram.png')
    return dend['ivl']

def heatmap_cluster(data_matrix,timepoint):
    """Produces a heatmap of the clustered count data"""
    D = np.zeros([np.shape(data_matrix)[0],np.shape(data_matrix)[0]])
    for i in range(np.shape(data_matrix)[0]):
        for j in range(np.shape(data_matrix)[0]):
            D[i,j] = abs(data_matrix[i] - data_matrix[j])**2  #use Wards method (other methods could be implemented here)
    fig = plt.figure()
    axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
    linked = sch.linkage(D, method='centroid')
    dend = sch.dendrogram(linked, orientation='right') # sets the oirentation root at the right
    axdendro.set_xticks([])
    axdendro.set_yticks([])
    #plot distance matrix
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
    index = dend['leaves']
    D=D[index,:]
    D=D[:,index]
    im = axmatrix.matshow(D, aspect='auto', origin='lower')
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    #plot color bar
    axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
    fig.colorbar(im, cax=axcolor)
    #display the heatmap
    fig.savefig(timepoint+'heatmap.png')

