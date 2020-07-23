## Tiedie result dir is $1

# input gene list file $1
databasepath="/export/home/kongyan/software/cancer_hall_mark_enrichment"


while read NAME;do cat $databasepath/database/hallmakrpathway/$NAME.txt $1|sort |uniq -d |xargs -n100|sed 's/ /,/g'; done < $databasepath/database/hallmakrpathway/name >tmp1
while read NAME;do cat $databasepath/database/hallmakrpathway/$NAME.txt $1|sort |uniq -d |wc -l; done < $databasepath/database/hallmakrpathway/name >tmp2
number=$(cat $1 $databasepath/database/hallmakrpathway/allgene |sort |uniq -d |wc -l)
awk '{print "'$number'"}' tmp1 >tmp3
awk '{print "3039"}' tmp1 >tmp4
paste  tmp2 tmp3 $databasepath/database/hallmakrpathway/number tmp4 tmp1 >tmp5
echo "data=read.table(file=\"tmp5\",header=F)
m=c()
for(i in 1:nrow(data))
{
m[i]=1-phyper(data[i,1],data[i,2],data[i,4],data[i,3],lower.tail = T)
}
write.table(file=\"0\",m)
" >tmp.R 
Rscript tmp.R
sed 1d 0|sed 's/\"//g;s/ /\t/g'|cut -f 2 >tmp6

paste $databasepath/database/hallmakrpathway/name $databasepath/database/hallmakrpathway/name2 tmp5 tmp6 |awk 'BEGIN{FS="\t"}{if($3>0){print $0"\t"($3/$4)/($5/$6)}}'|sed 's/\t/#/g'>main_enrich.xls
sed -i '1i Main Cancer hallmark category#Sub Hallmark pathway#my proteins in this hallmark#my proteins in all hallmark#genes in this hallmark#genes in all hallmark#gene list#p.value#enrichment ratio' main_enrich.xls
#rm tmp*
## Tiedie result dir is $1

#echo "Upper gene number is: "$a"; Downstream gene number is: "$b"; Linker gene number is: "$c
while read NAME;do cat $databasepath/database/hallmakrpathway/hallmark/$NAME $1|sort |uniq -d |xargs -n100|sed 's/ /,/g'; done < $databasepath/database/hallmakrpathway/hallmark/file >tmp1
while read NAME;do cat $databasepath/database/hallmakrpathway/hallmark/$NAME $1|sort |uniq -d |wc -l; done < $databasepath/database/hallmakrpathway/hallmark/file >tmp2
number=$(cat $1 $databasepath/database/hallmakrpathway/allgene |sort |uniq -d |wc -l)
awk '{print "'$number'"}' tmp1 >tmp3
awk '{print "3039"}' tmp1 >tmp4
paste  tmp2 tmp3 $databasepath/database/hallmakrpathway/hallmark/number tmp4 tmp1 |awk 'BEGIN{FS="\t"}{if(!$5){print $0"null";next}print $0}' >tmp5
Rscript tmp.R
sed 1d 0|sed 's/\"//g;s/ /\t/g'|cut -f 2 >tmp6

paste $databasepath/database/hallmakrpathway/hallmark/name $databasepath/database/hallmakrpathway/hallmark/name2 tmp5 tmp6 |awk 'BEGIN{FS="\t"}{print $0"\t"($3/$4)/($5/$6)}'|sed 's/\t/#/g'|cut -d "#" -f 1,2,4- >small_enrich.xls
sed -i '1i Main Cancer hall mark pathway#hallmark pathway#my proteins in this hallmark#my proteins in all hallmark#genes in this hallmark#genes in all hallmark#gene list#p.value#enrichment ratio' small_enrich.xls

#rm tmp*
echo "Name p1 p2" |sed 's/ /\t/g' >data.txt 
sed 's/#/\t/g' small_enrich.xls|cut -f 3,9  |sed 1d |awk 'BEGIN{FS="\t"}{print $1"\t0.001\t"(-1)*log($2)/log(10)}' >>data.txt
echo "data=read.table(file=\"data.txt\",header=T,sep=\"\t\")

library(reshape)
library(ggplot2)
library(plyr)
library(grDevices)

nba=data 
jpeg(\"wheel.jpg\",width=1200,height=1200) 

nba.m <- melt(nba)
nba.m <- ddply(nba.m, .(variable), transform, value = value)
nba.m\$var2 = as.numeric(nba.m\$variable) + 35
y_labels = levels(nba.m\$variable)
y_breaks = seq_along(y_labels) + 35

nba.labs <- subset(nba.m, variable==levels(nba.m\$variable)[nlevels(nba.m\$variable)])
nba.labs <- nba.labs[order(nba.labs\$Name),]
nba.labs\$ang <- seq(from=(360/nrow(nba.labs))/1.5, to=(1.5*(360/nrow(nba.labs)))-360, length.out=nrow(nba.labs))+80
nba.labs\$hjust <- 0
nba.labs\$hjust[which(nba.labs\$ang < -90)] <- 1
nba.labs\$ang[which(nba.labs\$ang < -90)] <- (180+nba.labs\$ang)[which(nba.labs\$ang < -90)]
p2 = ggplot(nba.m, aes(x=Name, y=var2, fill=value)) +
     geom_tile(colour=\"white\") +
     geom_text(data=nba.labs, aes(x=Name, y=var2+0.6,
        label=Name, angle=ang, hjust=hjust), size=4) +
     scale_fill_gradient(low = \"white\", high = \"darkred\") +
     ylim(c(0, max(nba.m\$var2) + 6.6)) +
     scale_y_discrete(breaks=y_breaks, labels=y_labels) +
     coord_polar(theta=\"x\") +
     theme(panel.background=element_blank(),
           axis.title=element_blank(),
           panel.grid=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks=element_blank(),
           axis.text.y=element_text(size=3))
print(p2)

dev.off() " >codewheel.R
xvfb-run Rscript codewheel.R

python merge.wheel.plot.py
rm tmp* 0
sed -i 's/#/\t/g' *.xls


mv cancer_hallmark_enrichment.jpg *.xls result/
rm *.jpg
rm data.txt
