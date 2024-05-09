#!/bin/bash 



region="";
# pop1_start=5;pop1_end=5;pop2_start=6;pop2_end=6;
pop1_start=5;pop1_end=9;pop2_start=10;pop2_end=14;
outfile=""

while getopts "r:a:b:c:d:o:" opt
do
  case $opt in
    r)
      region="-R $OPTARG"
      ;;
    a)
      pop1_start="$OPTARG"
      ;;
    b)
     pop1_end="$OPTARG"
      ;;
    c)
     pop2_start="$OPTARG"
      ;;
    d)
     pop2_end="$OPTARG"
      ;;     
    o)
     outfile="$OPTARG"
      ;;
    ?)
      echo "getopts param error"
      exit 1;;
  esac
done
shift $((OPTIND -1))

if [ $# -ne 3 ] ; then
echo "$0 vcffile [-a 5 -b 9 -c 10 -d 14 -r region_file -o pos.txt] sample1 sample2"
exit 0
fi

vcf=$1
sample1=$2
sample2=$3


read h1 h2 h3 h4 h5 h6 < <(bcftools view  -m 2 -M 2 $vcf  -s ${sample1},${sample2} $region|bcftools view - -g ^miss -v snps|bcftools query -f  '%CHROM  %POS  %REF  %ALT{0} [ %GT]\n'|grep -v '0/1'|awk \
-v pop1_start=$pop1_start  -v pop1_end=$pop1_end -v pop2_start=$pop2_start -v pop2_end=$pop2_end  -v outfile="$outfile" ' 
{ 
    cnt00=cnt11=0;for(i=pop1_start;i<=pop1_end;++i) {if($i=="0/0") cnt00++;else if($i=="1/1") cnt11++;  } if(cnt00==0 && cnt11==0) next; key=(cnt00 && cnt11)?"0/1":(cnt00?"0/0":"1/1");
    cnt00=cnt11=0;for(i=pop2_start;i<=pop2_end;++i) {if($i=="0/0") cnt00++;else if($i=="1/1") cnt11++;  } if(cnt00==0 && cnt11==0) next; key=sprintf("%s_%s",key,(cnt00 && cnt11)?"0/1":(cnt00?"0/0":"1/1")); 
    if(outfile) print $1,$2,substr(key,1,3),substr(key,5,3) > outfile;
    map[key]++; 
} 
END{
    close(outfile);
    h1=h2=h3=h4=0;
    for(x in map) {
          if(x=="0/0_0/0" || x=="1/1_1/1") continue; 
          if(x=="0/1_0/1") h4+=map[x]; 
          else if(x=="0/0_1/1" || x=="1/1_0/0") h3+=map[x]; 
          if(substr(x,1,3)=="0/1") h1+=map[x];
          if(substr(x,5,3)=="0/1") h2+=map[x];
    }
    print h1,h2,h3,h4,h1-h4,h2-h4;
}')

echo $sample1 $sample2 $h1 $h2 $h3 $h4 $h5 $h6
