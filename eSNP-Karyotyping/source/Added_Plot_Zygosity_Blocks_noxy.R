# Plot_Zygosity_Blocks_noxy
#
# Plot blocks of heterozygous and homozygous SNPs
# @param Table The deletion table containing the output of the DeletionTable function
# @param window The block size in bp, usually 1500000
# @param Max How many Heterozygous SNPs need to be in a block to get the full color, usually 6
# @param Max2 How many Homozygous SNPs need to be in a block to get the full color, usually 60
# @param Organism "Human" or "Mouse"
# @export
# @return None

Plot_Zygosity_Blocks_noxy<-function(Table,Window,Max,Max2,Organism){

  max=Max
  max2=Max2
  window=Window
  tbl=Table
  
  if(Organism == "Human"){
centromere_pos <- c(123400000,  93900000,  90900000,  50000000,  48800000,  59800000,  60100000,  45200000,
                    43000000,  61000000,  10400000,  39800000,  53400000,  35500000,  17700000,  17200000,
                    19000000,  36800000,  25100000,  18500000,  26200000,  28100000
                    #,  12000000,  15000000
                   )

chr_size <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
              133797422, 135086622, 133275309, 114364328, 107043718, 101991189,  90338345,  83257441,  80373285,
              58617616,  64444167,  46709983,  50818468
              #, 156040895,  57227415)
              )
    lb=c(seq(from = 150, to = 10, by = -20),seq(from = 10, to = 150, by = 20))
    plot(c(1:22),rep(0,22),ylim=c(-150000000,150000000),cex.axis=0.7,xlab = "Chromsome",ylab="Position (Mbp)",xaxt = "n",yaxt="n")
    axis(2, seq(from = 150000000, to = -150000000, by = -20000000), labels=lb,cex.axis=0.7,las=1)
    mx=22
  }
  
  if(Organism == "Mouse"){
    chr_size=c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768  ,94987271  ,90702639  ,61431566  ,171031299,91744698)
    centromere_pos=rep(200000000,21)
    lb=c(seq(from = 0, to = 200, by = 20))
    plot(c(1:21),rep(200000000,21),ylim=c(0,215000000),cex.axis=0.7,xlab = "Chromsome",ylab="Position (Mbp)",xaxt = "n",yaxt="n")
    axis(2, seq(from = 200000000, to = 0, by = -20000000), labels=lb,cex.axis=0.7,las=1)
    mx=21
  }
  
  tbl$snp=tbl$AF1>0
  
  
  total_homo=sum(tbl$snp==FALSE)
  total_hetro=sum(tbl$snp==TRUE)
  pval=NULL
  rat=NULL
  
  for (i in 1:mx){
    if (i<(mx-1)){
      chr=i
      tbl2=tbl[tbl$chr.x==chr,]
      
      d=tbl2[tbl2$start<centromere_pos[i],]
      homo=sum(d$snp==FALSE)
      hetro=sum(d$snp==TRUE)
      rat=c(rat,homo/hetro)
      
      d=tbl2[tbl2$start>centromere_pos[i],]
      homo=sum(d$snp==FALSE)
      hetro=sum(d$snp==TRUE)
      rat=c(rat,homo/hetro)
      
    }
  }
  
  
  rat[is.infinite(rat)==TRUE]=total_homo/total_hetro*5.5
  pval=NULL
  
  for (i in 1:((mx-2)*2)){
    np=1
    if(is.na(rat[i])==FALSE){np=t.test(rat,mu=rat[i])$p.value}
    pval = cbind(pval,np)
  }
  
  
  pval=p.adjust(pval,method  ="fdr")
  pval=log(pval,10)*(-1)
  rat[is.infinite(rat)==TRUE]=total_homo/total_hetro
  rat[is.na(rat)==TRUE]=total_homo/total_hetro
  
  y=pval>3 & rat>(total_homo/total_hetro*5)
  loc=1
  for (i in 1:(mx-2)){
    for (j in 1:2){
      
      if (y[loc]==FALSE){col="white"}
      if (y[loc]==TRUE){col=rgb(249,244,93,maxColorValue = 255)}
      if(j==1){rect(xleft = i-0.45,ybottom = +centromere_pos[i],xright = i+0.45,0,col = col,border = NA)}
      if(j==2){rect(xleft = i-0.45,ybottom = +centromere_pos[i]-chr_size[i],xright = i+0.45,0,col =col,border=NA)}
      loc=loc+1
      
    }
  }
  
  
  
  
  for (i in 1:mx){
    
    chr=i
    tbl2=tbl[tbl$chr.x==chr,]
    tr=tbl2[tbl2$snp==T,]
    fl=tbl2[tbl2$snp==F,]
    
    
    win=ceiling(chr_size[i]/window)
    pos=centromere_pos[i]
    library(gplots)
    cl=colorpanel(max, "grey", "red")
    cl2=colorpanel(max2, "grey", "blue")
    
    for (j in 1:win){
      top=j*window
      bottom=(j-1)*window
      num_fl=sum(fl$start>=bottom & fl$start<=top)
      num_tr=sum(tr$start>=bottom & tr$start<=top)
      
      if (num_tr>max){num_tr=max}
      if (num_fl>max2){num_fl=max2}
      
      rect(i,pos,i+0.4,pos-window,col=cl[num_tr],border = NA)
      rect(i,pos,i-0.4,pos-window,col=cl2[num_fl],border = NA)
      
      
      pos=pos-window
    }
    
    
    
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=10,col="black")  
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=8,col="gray50")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=7,col="gray53")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=6,col="gray59")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=4,col="gray75")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=2,col="gray85")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=1,col="gray90")
    if(Organism=="Human"){
      points(i,0,pch=16,col="grey13")
      if (i<23) {text(i,centromere_pos[i]+18000000,i,cex=0.8)}
      #if (i==23) {text(i,centromere_pos[i]+18000000,"X",cex=0.8)}
      #if (i==24) {text(i,centromere_pos[i]+18000000,"Y",cex=0.8)}
    }
    
    if(Organism=="Mouse"){
      points(i,200000000,pch=16,col="grey13")
      if (i<20) {text(i,212000000,i,cex=0.8)}
      if (i==20) {text(i,212000000,"X",cex=0.8)}
      if (i==21) {text(i,212000000,"Y",cex=0.8)}}
    
  }
  
  
}