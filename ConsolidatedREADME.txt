Consolidated analyses for Patellifolia poolseq study

This document offers details on the analyses conducted, but scripts are not intended to be
used as is to perform analyses in another system.  System dependencies would necessitate
many changes to the code here in order to adapt the pipeline.

In what follows the integers 50-55 are used to represent the Patellifolia pools. 50-52 are
P. patellaris, 53-54 P. procumbens, 55 P. webbiana. Analyses were conducted on two
different HPC resources, named 'ceres' and 'blip' so some paths or variables may refer to
those by name.  Analyses were submitted to the HPC using SLURM, or conducted in interactive
sessions, often using GNU Parallel.

 



### GENOME ASSEMBLY ###
#Notes on de novo assembly of 6 Patellifolia genomes using Masurca

#De novo genome assembly was performed on a single compute node containing two 10 core
#Intel Xeon E5-2680 v2 @ 2.80GHz processors and 512 GB RAM

#Edit the default config file called sr_config_example.txt with the following changes
#relative to default:


PE= pe 150 20 forward.fastq.gz reverse.fastq.gz
#JUMP= sh 3600 200 /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
#PACBIO=/FULL_PATH/pacbio.fa
#OTHER=/FULL_PATH/file.frg
USE_LINKING_MATES = 1
NUM_THREADS = 40
JF_SIZE = 5000000000


#Run the following to produce the assemble.sh file (specify full paths):

/home/reevesp/bin/MaSuRCA-3.2.4/bin/masurca /home/reevesp/bin/MaSuRCA-3.2.4/50/config50.txt; #this command saves assemble.sh to the current directory

#In assemble.sh, change the line controlling super reads, adding '-low-memory', as in the 
#following diff output:
#This will make it cache to disk during the 'Computing super reads from PE' step.


< createSuperReadsForDirectory.perl -low-memory -l $KMER -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.pe.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.fasta -t 40 -mikedebug work1 pe.cor.fa 1> super1.err 2>&1
---
> createSuperReadsForDirectory.perl -l $KMER -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.pe.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.fasta -t 40 -mikedebug work1 pe.cor.fa 1> super1.err 2>&1

sed -i 's/createSuperReadsForDirectory.perl/createSuperReadsForDirectory.perl -low-memory/' assemble.sh;


#In assemble.sh, enable parallelism by uncommenting the following lines like:
# To run tasks in parallel
run_bg () {
  semaphore -j $NUM_THREADS --id masurca_$$ -- "$@"
}
run_wait () {
  semaphore -j $NUM_THREADS --id masurca_$$ --wait
}

#Run assemble.sh
./assemble.sh 2>&1 | tee outputlog.txt; #write stderr and stdout to a log that can be monitored, useful when running in screen





### INDEX GENOME ASSEMBLY ###
#crude genome assemblies, hereafter also referred to as "reference", from masurca were
#named 5[0-5]Hs1pro1REVdref.fasta for the 6 pools
for i in {50..55};
  do echo "$i";
    samtools faidx "$i"Hs1pro1REVdref.fasta; #fai file contains contig names and lengths
 done;
seq 50 1 55 | parallel 'bwa index '{}'Hs1pro1REVdref.fasta';





### TRIM RAW SEQUENCE READS ###
mypp() {
        i=$1;
        echo "$i";
        java -jar ~/bin/trimmomatic.jar PE -threads 40 17134D-01-$i*_L001_R1_001.fastq.gz 17134D-01-$i*_L001_R2_001.fastq.gz \
          "$i"output_forward_paired.fq.gz "$i"output_forward_unpaired.fq.gz "$i"output_reverse_paired.fq.gz "$i"output_reverse_unpaired.fq.gz \
          ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50;
}
export -f mypp;
seq 50 55 | parallel --env mypp mypp; #use parallel because multithreading doesn't work.





### MAP READS TO POOL-SPECIFIC REFERENCE ###
#Align reads from each pool to its reference assembly using bwa-mem, and markdups using sambamba
#This saves all output, intermediate and otherwise, to the local node. Things can get messy if it fails.
mypp() { 
       i=$1;
       f="$wd"/"$i"output_forward_paired.fq.gz; #forward reads
       r="$wd"/"$i"output_reverse_paired.fq.gz; #reverse reads
       v=$(echo "$wd" | sed 's/data/ref/')/"$i"Hs1pro1REVdref.fasta; #path to reference assembly
       thr=$(lscpu | grep "^CPU(s):" | awk '{print $2}'); #max threads
       tmpd="/state/partition1/tmpbwa$i";
       if [ ! -d "$tmpd" ]; then mkdir "$tmpd"; fi; #make a local tmp directory on the main disk
       readgroup="@RG\tID:$i\tSM:$i\tPL:illumina\tLB:na\tPU:na";

       #do bwa-mem alignment, save result to local directory
       #stdout from bwa mem goes to samtools sort, a coordinate sorted bam file is the result
       # -k 19 is the default seed length, it can be varied below
       /home/reevesp/bin/bwa mem -R "$readgroup" -t "$thr" -k 19 "$v" "$f" "$r" | /home/reevesp/bin/samtools sort -O BAM --threads "$thr" -T "$tmpd"/tmp.$i.bam -o "$tmpd"/$i"_aligned_reads.Hs1pro1.bam"; 

       #markdups and clean up
       /home/reevesp/bin/sambamba markdup -t "$thr" --tmpdir="$tmpd" --overflow-list-size 6000000 --remove-duplicates "$tmpd"/$i"_aligned_reads.Hs1pro1.bam" "$tmpd"/$i"_markdups.Hs1pro1.bam";
       rm "$tmpd"/$i"_aligned_reads.Hs1pro1.bam";
}
export -f mypp;
cd /share/space/reevesp/patellifolia/data;
wd=$(pwd); export wd;
seq 50 1 55 | parallel --sshloginfile ~/machines --env wd --env mypp mypp;





### PHASE READS ###
#The below code can be pasted to a file hapxsubmitter.slm, which is submitted as a slurm job.
#Because hundreds of thousands of jobs will be submitted by hapxsubmitter.slm (one per
#contig in the crude reference assembly), much code is devoted to controlling that process, including
#keeping the squeue of a manageable length for other users, and consolidating output files
#so as not to overwhelm the filesystem.

#The script calls hapx (https://github.com/NCGRP/hapx), which calls bwa-mem and samtools and makes
#extensive use of GNU parallel.  All of those software are dependencies.

#You need to edit the #SBATCH parameters in the enveloping slurm batch file hapxsubmitter.slm.
#You also need to edit the $clus value and others for running on blip or ceres prior to submission.


#####BEGIN SLURM#####
#!/bin/bash
#SBATCH --job-name=hapxsub #name of the job submitted
#SBATCH -p long60 #name of the queue you are submitting to, also scavenger,scavenger-mem,scavenger-mem768
#SBATCH -N 1 #number of nodes in this job
#SBATCH -n 1 #number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH --mem=2G
#SBATCH -t 60-00:00:00 #time allocated for this job hours:mins:seconds
#SBATCH -o "/dev/null" # standard out goes nowhere
#SBATCH -e "stderr.hapxsub.%j.%N.%A.%a" #optional but it prints out standard error

module load samtools;

#modify these values depending on HPC you are using
clus=ceres; #define cluster to create slurm file for, "ceres" "blip"
memz=3780M; #memory for slurm, default this setting to "6G" for blip, "3780M" for ceres
userz="reev" #"ann.reil"  "reev" #this is a truncation of the user name, grepable from squeue
nodez="scavenger"; #"scavenger" "medium,long,mem,mem768" "CLUSTER"
poolz="51" #a single pool name, or a string like: "50 51 52", for use in a standard bash for loop to iterate over pools

if [[ "$clus" == "blip" ]];
then ppath="/share/space/reevesp/patellifolia/";
  nl=1000;  #1000, number of jobs to submit to slurm before waiting for them to finish/testing whether enough have finished to merge
  nr=250; #250, number of jobs in slurm queue, below which more jobs are submitted
  ml=800; #800, merge level, number of completed jobs before a merge will be performed, ulimit is 1024 which will sometimes cause samtools to fail if this is set to 1000
elif [[ "$clus" == "ceres" ]];
then ppath="/lustre/project/patellifolia/";
  nl=50; #20
  nr=250; #189
  ml=1000; #1000
fi;
    
for jk in $poolz;
  do >runlog"$jk".txt;
    
    rref="$ppath""ref/"$jk"Hs1pro1REVdref.fasta";
    bbam="$ppath""map/"$jk"_markdups.Hs1pro1.bam"; #alignment file
    
    #set up to pass $nl lines at a time to slurm
    st=$(seq 1 "$nl" $(wc -l "$rref".fai | cut -d' ' -f1));
    en=$(seq "$nl" "$nl" $(wc -l "$rref".fai | cut -d' ' -f1));
    en="$en"$'\n'$(wc -l "$rref".fai | cut -d' ' -f1); #add final line of fai file
    ranges=$(paste -d':' <(echo "$st") <(echo "$en") | tr '\n' ' ');
    
    #reverse order of lines in $ranges so longest contigs run first
    #this was for testing purposes and doesn't really do what is advertised
    r1=$(echo "$ranges" | tr " " "\n" | awk 'NF' | tac | tr "\n" " ");
    ranges="$r1";
    
    jobticker=1; #counts the number of batches submitted
    for k in $ranges;
      do s=$(echo "$k" | cut -d: -f1);
        e=$(echo "$k" | cut -d: -f2);
        a=$(sed -n "$s"','"$e"'p;'"$e"'q' "$rref".fai | cut -d$'\t' -f1,2 | sed 's/\t/:1-/'); #assemble the $nl line set of contig:site-ranges
    
      #submit to slurm
      for aa in $a;
        do echo "submitting $aa range $k";
          { echo -n "submitting $aa range $k ";
          outfol=$(sed 's/:/_/' <<<"$aa"); #generate the output folder name
          mkdir "$jk""$outfol"; #set up the output folder in advance so you can write the stderr to it
          
          #create slurm file
          if [[ "$clus" == "blip" ]];
          then b='#!/bin/bash
#SBATCH --job-name='"$aa"' #name of the job submitted
#SBATCH -p '"$nodez"' #short,medium,long,mem,mem768 #name of the queue you are submitting job to
#SBATCH -N 1 #number of nodes in this job
#SBATCH -n 1 #number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH --mem='"$memz"'
#SBATCH -t 7-00:00:00 #time allocated for this job hours:mins:seconds
#SBATCH -o '"/dev/null"' # standard out goes nowhere
#SBATCH -e '"$jk$outfol/stderr.%j.%N"' #optional but it prints out standard error

# submit this array like: sbatch <(echo "$b")

#execute hapx
hapx.sh -r '"$rref"' \
    -b '"$bbam"' \
    -o '"$jk$outfol"' -mb -q 1 -sp -s <(echo '"$aa"');
    
#log job completion
echo '"completing $aa range $k Completed batch job "'${SLURM_JOB_ID} >> '"$jk$outfol"'/stderr.${SLURM_JOB_ID}.${SLURMD_NODENAME};

#End of file
';
 
           else b='#!/bin/bash
#SBATCH --job-name='"$aa"' #name of the job submitted
#SBATCH -p '"$nodez"' #scavenger #medium,long,mem,mem768 #name of the queue you are submitting to, also scavenger,scavenger-mem,scavenger-mem768
#SBATCH -N 1 #number of nodes in this job
#SBATCH -n 1 #number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH --mem='"$memz"'
#SBATCH -t 7-00:00:00 #time allocated for this job hours:mins:seconds
#SBATCH -o '"/dev/null"' # standard out goes nowhere
#SBATCH -e '"$jk$outfol/stderr.%j.%N"' #optional but it prints out standard error

# submit this array like: sbatch <(echo "$b")

# add regular job commands like module load and running scientific software
module load parallel/20151022;
module load bwa;
module load samtools;
module load muscle;
#module load gcc/5.3.0

#execute hapx
/project/patellifolia/bin/hapx.sh -r '"$rref"' \
    -b '"$bbam"' \
    -o '"$jk$outfol"' -mb -q 1 -sp -s <(echo '"$aa"');

#log job completion
echo '"completing $aa range $k Completed batch job "'${SLURM_JOB_ID} >> '"$jk$outfol"'/stderr.${SLURM_JOB_ID}.${SLURMD_NODENAME};

#End of file
';
         fi; #if clus="blip" or "ceres"
    
          #submit to slurm queue
          sbatch <(echo "$b"); #use process substitution to submit variable $b as if it is a slurm batch file
          } >> runlog"$jk".txt; #write progress to log
  
        done; #for aa in $a
        
        #Test whether enough jobs have completed to perform a merge, if not, don't merge
        #and test whether queue is full. If so, merge sorted bam files containing flashed
        #read pairs for each contig, also merge reads in global.fa files, then test whether
        #queue is full       
        sleep 10; #wait while slurm queues up
        while true;
        do af=$(find */alignments -name "*.bam"); #list of bam files at current time
          bamcount=$(echo "$af" | wc -l);
          echo "$bamcount bam files present." >> runlog"$jk".txt;
          if (( "$bamcount" > "$ml" ));
          then
            { echo "  Merging job $jobticker.";
            samtools merge -b <(echo "$af") "$jobticker".bam;
            echo "$af" | sed 's/_aligned_haps.bam/.global.fa/' | xargs cat > "$jobticker".fa;
        
            #clean up
            echo "$af" | cut -d'/' -f1 | xargs -i -P24 sh -c 'rm -r *{}*'; #delete all folders from the current processed set $af
            jobticker=$(( $jobticker + 1 ));
            } >> runlog"$jk".txt;
          fi;
          #keep queue full, test whether new jobs need to be submitted
          if [[ "$nodez" == "scavenger" ]];
          then qn=$(squeue | awk -v userz="$userz" '$4~userz{print $0}' | grep scavenger | grep -v hapxsub | grep -v " CG " | wc -l); #number of jobs in queue scavenger queue
          else qn=$(squeue | awk -v userz="$userz" '$4~userz{print $0}' | grep -v scavenger | grep -v hapxsub | grep -v " CG " | wc -l); #number of jobs in non-scavenger queue
          fi;
          
          if (( "$qn" > "$nr" ));
            then echo "Running $s-$e, $qn jobs in queue          "$(date) >> runlog"$jk".txt;
              sleep 10;
            else break; #break allows more jobs to be submitted
          fi;
        done;
        
      done; #k in $ranges
      


      #final processing step for lingering runs after all $ranges have been processed
      while true;
      do if [[ "$nodez" == "scavenger" ]];
          then qn=$(squeue | awk -v userz="$userz" '$4~userz{print $0}' | grep scavenger | grep -v hapxsub | grep -v " CG " | wc -l); #number of jobs in queue scavenger queue
          else qn=$(squeue | awk -v userz="$userz" '$4~userz{print $0}' | grep -v scavenger | grep -v hapxsub | grep -v " CG " | wc -l); #number of jobs in non-scavenger queue
        fi;
      if (( "$qn" > 1 ));
         then echo "Running $s-$e, $qn jobs remaining          "$(date) >> runlog"$jk".txt;
           sleep 60;
         else echo "Final processing." >> runlog"$jk".txt;
           break;
         fi;
      done;

      #merge sorted bam files containing flashed read pairs for each contig, also merge reads in global.fa files
      { echo "Merging job $jobticker.";
      af=$(find */alignments -name "*.bam"); #all bam files in output folder at this moment (changes constantly due to slurm)
      samtools merge -b <(echo "$af") "$jobticker".bam;
      echo "$af" | sed 's/_aligned_haps.bam/.global.fa/' | xargs cat > "$jobticker".fa;
      
      #clean up
      echo "$af" | cut -d'/' -f1 | xargs -i -P24 sh -c 'rm -r *{}*'; #delete all folders from the current processed set $af
      jobticker=$(( $jobticker + 1 ));
      } >> runlog"$jk".txt;

  done; #for i in list of pool numbers

#####END SLURM#####





### POST RUN VALIDATION ###

## FOR HPC BLIP (HPC CERES BELOW) ##

#Consolidate all results for the run in a folder named like 55fra and move to folder
#/share/space/reevesp/patellifolia/FRAresults.
#Values 5[0-5] will be observed in the following code and represent codes for each of the 6 pools

#check for runs where no reads aligned, move them to a separate folder
cd /share/space/reevesp/patellifolia/FRAresults/55fra;
mkdir NoReads;
find . -name "NoReadsSoNoAlignmentPossible" | cut -d'/' -f1-2 | xargs mv -t NoReads;

#identify the types of errors that have occurred in the remaining run folders
b=$(ls -d 50jcf*); #get list of remaining runs
for bb in $b;
  do echo "$bb";
    grep -v "0 sec " "$bb"/stderr* | grep -v '/dev/tty'; #show lines that are not progress bar
    grep " DUE TO TIME LIMIT" "$bb"/stderr*;
    echo;
  done;

#if some jobs exceeded memory limit, set memory from 8GB to 48GB and try all jobs again on compute-0-9 only
cd /share/space/reevesp/patellifolia/FRAresults/55fra;
a=$(ls -d 55jcf* | sed 's/_/:/g' | sed 's/^55//'); #get list of remaining runs from e.g. /share/space/reevesp/patellifolia/FRAresults/55fra

cd /share/space/reevesp/patellifolia/flashedreadarchive; #this folder must be empty when you start
su; #use root to turn off all nodes but compute-0-9
scontrol update nodename=compute-0-[0-5,7,10] state=drain reason="rerun55";
exit;

#use the following steps to perform the resubmission of failed jobs on blip
clus=blip; #define cluster to create slurm file for, "ceres" "blip"
memz=48G; #memory for slurm, default this setting to "8G" for blip, "3780M" for ceres
ppath="/share/space/reevesp/patellifolia/";
nl=$(echo "$a" | wc -l);  #number of jobs to submit to slurm before waiting for them to finish/testing whether enough have finished to merge
nr=232; #number of jobs in slurm queue, below which more jobs are submitted
ml=$(echo "$a" | wc -l); #merge level, number of completed jobs before a merge will be performed

jk=55;
rref="$ppath""ref/"$jk"Hs1pro1REVdref.fasta";
bbam="$ppath""map/"$jk"_markdups.Hs1pro1.bam"; #alignment file
k="rerun";

>runlog"$jk".txt;

#In the SLURM script for section ### PHASE READS ###, paste a portion starting at 'for aa in $a;' 
#and ending before 'done; #for aa in $a' to run loop to resubmit jobs to compute-0-9.
#This will perform the basic runs, then merge manually:
{ echo "  Merging job $jobticker.";
af=$(find */alignments -name "*.bam"); #all bam files in output folder at this moment (changes constantly due to slurm)
samtools merge -b <(echo "$af") "$jobticker".bam;
echo "$af" | sed 's/_aligned_haps.bam/.global.fa/' | xargs cat > "$jobticker".fa;

#clean up
echo "$af" | cut -d'/' -f1 | xargs -i -P24 sh -c 'rm -r *{}*'; #delete all folders from the current processed set $af
jobticker=$(( $jobticker + 1 ));
} >> runlog"$jk".txt;


#archive contigs that did not finish or had errors twice within the 7 day max time allotted
cd /share/space/reevesp/patellifolia/FRAresults/55fra;
mkdir FailedContigs;
cd /share/space/reevesp/patellifolia/flashedreadarchive;
for i in $(ls -d 55*);
  do echo $i;
    a=$(find $i/alignments -type f | wc -l); #find alignments folders with 0 files within
    if (( $a == 0 ));
    then mv $i /share/space/reevesp/patellifolia/FRAresults/55fra/FailedContigs;
    fi;
  done;

#consolidate the newly redone steps with the files created during the first run
find . -name "*.fa" -exec cp {} /share/space/reevesp/patellifolia/FRAresults/55fra \;
find . -name "*.bam" -exec cp {} /share/space/reevesp/patellifolia/FRAresults/55fra \;

#combine *.fa and *.bam files into 1
cd /share/space/reevesp/patellifolia/FRAresults/55fra;
mkdir 55fraFinal;
cp /share/space/reevesp/patellifolia/ref/55Hs1pro1REVdref.fasta 55fraFinal;
cat *.fa > 55fraFinal/55fra.fa;
sambamba merge -t24 55fraFinal/55fra.bam *.bam;
samtools flagstat 55fraFinal/55fra.bam > 55fraFinal/55fra.bam.flagstat; #show some stats about mapped reads

#make sure all indexes are created for references
cd /share/space/reevesp/patellifolia/ref;
mypp() {
       i=$1;
       echo "$i";
       makeblastdb -in "$i"Hs1pro1REVdref.fasta -parse_seqids -dbtype nucl;
       #samtools faidx "$i"Hs1pro1REVdref.fasta;
       #bwa index "$i"Hs1pro1REVdref.fasta;
}
export -f mypp;
seq 50 1 55 | parallel --env mypp mypp;

#clean up directory on cluster, transfer to NAS for long term storage, compress FRA
cd /share/space/reevesp/patellifolia/FRAresults/55fra;
rm -r 55jcf*; #remove directories for runs that failed first round
mkdir ConstituentResults;
mv *.fa ConstituentResults;
mv *.bam ConstituentResults;
mv *runlog.txt ConstituentResults;
cd ../;

tar -zcvf 55fraFinal.tar.gz 55fraFinal; #compresses 77GB final fra to 27GB. raw reads = 30GB





### POST RUN VALIDATION ###

## FOR HPC CERES (HPC BLIP ABOVE) ##

#consolidate all results for the run in a folder named like 55fra and move to folder
#place all combined results files into a folder, later they will all be combined again

#check for runs where no reads aligned, move them to a separate folder
r=53; #50 51 52 53 54
cd ~/patellifolia/FRAresults;
mkdir "$r"fra; 
cd ~/patellifolia/FRAresults/"$r"fra;
mv ~/patellifolia/flashedreadarchive"$r"/* ~/patellifolia/FRAresults/"$r"fra; #move results
mkdir ConstituentResults;
mv *.bam ConstituentResults;
mv *.fa ConstituentResults;

mkdir NoReads;
find . -name "NoReadsSoNoAlignmentPossible" | cut -d'/' -f1-2 | xargs mv -t NoReads;

#check for runs with *.bam and *.fa files but which, for some reason, were not combined
#by the submission script
find "$r"*jcf*/alignments -name "*.bam" | wc -l; #confirm that both bam and fa files are present
find "$r"*jcf*/alignments -name "*.fa" | wc -l;
c=$(find "$r"*jcf*/alignments -name "*.bam" | cut -d'/' -f1); #list of completed but uncombined runs

#copy completed but uncombined result files to ConstituentResults
find "$r"*jcf*/alignments -name "*.bam" | xargs -I{} cp {} ConstituentResults;
find "$r"*jcf*/alignments -name "*.fa" | xargs -I{} cp {} ConstituentResults;

#archive runs with uncombined results
mkdir UncombinedResultsRuns;
echo "$c" | xargs -I{} mv {} UncombinedResultsRuns;

#identify the types of errors that have occurred in the remaining run folders
#in general, ceres jobs will not be memory limited, but may be time limited or just not run
b=$(ls -d "$r"*jcf*); #get list of remaining runs
for bb in $b;
  do echo "$bb";
    grep -v "0 sec " "$bb"/stderr* | grep -v '/dev/tty'; #show lines that are not progress bar
    grep " DUE TO TIME LIMIT" "$bb"/stderr*;
    echo;
  done;


#use the following steps to perform the resubmission of failed jobs on ceres
cd /home/pat.reeves/patellifolia/FRAresults/"$r"fra;

a=$(ls -d "$r"*jcf* | sed 's/_/:/g' | sed 's/^'$r'//'); #get list of remaining runs from e.g. /share/space/reevesp/patellifolia/FRAresults/55fra
cd /home/pat.reeves/patellifolia/flashedreadarchive"$r"; #this folder must be empty when you start

clus=ceres; #define cluster to create slurm file for, "ceres" "blip"
memz=3780M; #memory for slurm, default this setting to "8G" for blip, "3780M" for ceres
nodez="scavenger"; #"scavenger" "medium,long,mem,mem768" "CLUSTER"
ppath="/lustre/project/patellifolia/";
nl=$(echo "$a" | wc -l);  #number of jobs to submit to slurm before waiting for them to finish/testing whether enough have finished to merge
nr=250; #number of jobs in slurm queue, below which more jobs are submitted
ml=$(echo "$a" | wc -l); #merge level, number of completed jobs before a merge will be performed

jk="$r";
rref="$ppath""ref/"$jk"Hs1pro1REVdref.fasta";
bbam="$ppath""map/"$jk"_markdups.Hs1pro1.bam"; #alignment file
k="rerun";

>runlog"$jk".txt;

#In the SLURM script for section ### PHASE READS ###, paste a portion starting at 'for aa in $a;' 
#and ending before 'done; #for aa in $a' to run loop to resubmit jobs to compute-0-9.
#This will perform the basic runs, then merge manually:
sshort; #start up 40 core node on ceres
module load samtools;
r=53; #50 51 52 53 54
cd /home/pat.reeves/patellifolia/flashedreadarchive"$r"; #this is the folder where you've just repeated failed runs from first round
jobticker="$r""redo";

af=$(find */alignments -name "*.bam"); #all bam files in output folder at this moment (changes constantly due to slurm)
samtools merge -b <(echo "$af") "$jobticker".bam;
echo "$af" | sed 's/_aligned_haps.bam/.global.fa/' | xargs cat > "$jobticker".fa;

#clean up
echo "$af" | cut -d'/' -f1 | xargs -i -P40 sh -c 'rm -r *{}*'; #delete all folders from the current processed set $af

#archive contigs that did not finish or had errors twice within the 7 day max time allotted
cd /home/pat.reeves/patellifolia/FRAresults/"$r"fra;
mkdir FailedContigs;
cd /home/pat.reeves/patellifolia/flashedreadarchive"$r";
for i in $(ls -d "$r"*/);
  do echo $i;
    a=$(find "$i"alignments -type f | wc -l); #find alignments folders with 0 files within
    if (( $a == 0 ));
    then mv $i /home/pat.reeves/patellifolia/FRAresults/"$r"fra/FailedContigs;
    fi;
  done;

#consolidate the newly redone steps with the files created during the first run
find . -name "*.fa" -exec cp {} /home/pat.reeves/patellifolia/FRAresults/"$r"fra/ConstituentResults \;
find . -name "*.bam" -exec cp {} /home/pat.reeves/patellifolia/FRAresults/"$r"fra/ConstituentResults \;

#make sure all indexes are created for references
#only need to do this once, reference transfer step below relies on it having been already completed
module load blast+;
cd /home/pat.reeves/patellifolia/ref;
mypp() {
       i=$1;
       echo "$i";
       makeblastdb -in "$i"Hs1pro1REVdref.fasta -parse_seqids -dbtype nucl;
       #samtools faidx "$i"Hs1pro1REVdref.fasta;
       #bwa index "$i"Hs1pro1REVdref.fasta;
}
export -f mypp;
seq 50 1 54 | parallel --env mypp mypp;

#combine *.fa and *.bam files into 1
cd /home/pat.reeves/patellifolia/FRAresults/"$r"fra;
mkdir "$r"fraFinal;
cp /home/pat.reeves/patellifolia/ref/"$r"Hs1pro1REVdref.fasta* "$r"fraFinal;
cat ConstituentResults/*.fa > "$r"fraFinal/"$r"fra.fa;
sambamba merge -t40 "$r"fraFinal/"$r"fra.bam ConstituentResults/*.bam;
samtools flagstat "$r"fraFinal/"$r"fra.bam > "$r"fraFinal/"$r"fra.bam.flagstat; #show some stats about mapped reads

#clean up directory on cluster, transfer to NAS for long term storage
cd /home/pat.reeves/patellifolia/FRAresults/"$r"fra;
rm -r "$r"jcf*; #remove directories for runs that failed first round
mkdir RuntimeFiles;
mv runlog*.txt RuntimeFiles;
mv stderr* RuntimeFiles;
mv hapxsubmitter"$r".slm RuntimeFiles;





### SPECIFICATION AND FINALIZATION OF ARCHIVE ###

#A phased read archive should contain:
#1. reference sequence + index files including bwa index (bwt,amb,ann,pac,sa), samtools faidx (fai)
#2. merged bam file + index that aligns phased reads to reference
#3. sorted concatenated fasta file of phased reads that map onto reference
#4. blast formatted db (nhr,nin,nog,nsd,nsi,nsq,nal) of (1) and (3)

#Make single line fasta and sort *fra.fa files for use with sgrep downstream (~40 minutes per sort)
for i in 50 51 52 53 54 55;
  do echo "$i";
    (LC_ALL=C sed -e '/^>/s/$/@/' -e 's/^>/#>/' < "$i"fra.fa | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d' && echo) | LC_ALL=C sort -T /scratch/tmp > /scratch/reevesp/"$i"frasorted.fa; #make output 1 line per sequence
  done;

#Make blastdbs for fasta files of flashed reads, takes about 6 hours total
pd=$(pwd);
seq 50 1 55 | parallel 'makeblastdb -in '"$pd"'/{}frablastdb/{}frasorted.fa -parse_seqids -dbtype nucl;'

#backup FRAs onto NAS at /share/Public/Data/PatReeves/PatellifoliaIlluminaAnalysis/reeves/FlashedReadArchive




### CALCULATE QUALITY OF FRA ###
#calculate phased read count
for i in $(seq 50 1 55);
  do echo -n "$i ";
    grep ^'>'rp "$i"frasorted.fa | wc -l;
  done;
			poolID, number of phased reads in FRA
			50 242901801
			51 144535771
			52 184079155
			53 258249687
			54 173652412
			55 179541876

#calculate total phased read length
mypp() {
       i=$1;
       echo -n "$i ";
       grep ^'>'rp "$i"frasorted.fa | cut -d' ' -f2 | tr -d '\n' | wc -c;
}
export -f mypp;

seq 50 1 55 | parallel --env mypp mypp;

			poolID, cumulative length of phased reads
			50 66756924241
			51 39401459615
			52 50336787189
			53 69199316527
			54 47342488152
			55 45140337698




### IDENTIFICATION OF SUGAR BEET EL10 TRANSCRIPTOME (GENIC) HOMOLOGS IN PATELLIFOLIA ###

#perform blastn search in parallel on blip node compute-0-12 using 24255 EL10 primary transcripts as query with basic blast settings:
#Bvulgaris_548_EL10_1.0.cds_primaryTranscriptOnly.fa,-num_threads 9 -num_descriptions 100000 -num_alignments 0 (10 hrs for all 6 pools compute-0-12)
ssh compute-0-12;
cd /scratch/reevesp/patellifolia/blastdbflashedreads; #starting folder contains blast dbs for phased reads with local path like ./5[0-5]frablastdb/5[0-5]fra.fa as in -db line in blastn command below

g="Bvulgaris_548_EL10_1.0.cds_primaryTranscriptOnly.fa"; #query file
pqf="/scratch/reevesp/patellifolia/EL10ref"; #path to folder with query file
time seq 50 1 55 | parallel 'time blastn -num_threads 9 -num_descriptions 100000 -num_alignments 0 \
                                    -db {}frablastdb/{}fra.fa \
                                    -query '"$pqf/$g"' \
                                    -out {}frablastdb_'"$g"'.txt';


#parse blastn output into individual files for each query, and count reads mapped
mkdir $(seq 50 1 55);
j="frablastdb_Bvulgaris_548_EL10_1.0.cds_primaryTranscriptOnly.fa.txt"; #suffix for blast output file
export j;
for i in $(seq 50 1 55); do mv "$i""$j" "$i"; done;

mypp() {  
       i=$1;
       echo $i;
       cd "$pd"/"$i";
       csplit --quiet "$i""$j" '/Query/' '{*}'; #split on 'Query' to files named xx##
       time find . -name "xx*" -print0 | xargs -0 -I{} sh -c 'echo -n {}" "; \
                                                  echo -n $(head -1 {} | cut -d" " -f2)" "; \
                                                  grep ^[rj][pc] {} | wc -l' > "$i"linecount.txt; #~8min
}
export -f mypp

pd=$(pwd); export pd;
seq 50 1 55 | parallel --env mypp --env pd --env j mypp; #~6.5min


#try to figure out if there is some kind of natural breakpoint in hit count
ssh compute-0-12;
cd /scratch/reevesp/patellifolia/blastdbflashedreads;
cat 5*/5*linecount.txt | sort -t' ' -k2,2 -k3,3nr | grep -v '\./xx00 ' > readcounts.txt;
#plot hit x frequency to identify a reasonable range of read-hits (40-1000), see Excel file plotreadcounts.xlsx a.k.a. Supplemental Figure 1
cut -d' ' -f3 readcounts.txt | sort | uniq -c | sed 's/^ *//g' | sort -t' ' -k2,2nr | awk -F' ' '{print $2,$1}' | tr " " "\t" > plotrc.txt;

#get gene names where 40-1000 reads were matched
mkdir b100000yescp; #make a new directory to hold continuing analyses
mv 5[0-5] b100000yescp;
for i in $(seq 50 1 55);
  do awk -F' ' '$3>39{print $0}' b100000yescp/"$i"/"$i"linecount.txt | awk -F' ' '$3<1001{print $0}' > b100000yescp/"$i"/"$i"goodlines.txt;
  done;

wc -l b100000yescp/*/*goodlines.txt; #number of good genes
			  16779 b100000yescp/50/50goodlines.txt
			  16784 b100000yescp/51/51goodlines.txt
			  16995 b100000yescp/52/52goodlines.txt
			  16909 b100000yescp/53/53goodlines.txt
			  17104 b100000yescp/54/54goodlines.txt
			  17179 b100000yescp/55/55goodlines.txt

#move the files with a reasonable number of reads to a folder 'goodgenes'
cd /scratch/reevesp/patellifolia/blastdbflashedreads/b100000yescp;
for i in $(seq 50 1 55);
  do echo "$i";
    mkdir "$i"/"$i"goodgenes; #to hold files listing reads mapped to genes with >40 <1000 reads
    cut -d' ' -f1 "$i"/"$i"goodlines.txt | sed 's:^\./::g' | parallel 'mv '$i'/{} '$i'/'$i'goodgenes';
    mkdir "$i"/"$i"repetgenes; #to hold files listing reads mapped to genes with <40 >1000 reads
    mv "$i"/xx[0-9]* "$i"/"$i"repetgenes;
  done;


#What proportion of genes in EL10 map to unique contigs in patellifolia?
#These will be called 'single copy orthologs' and are determined for each pool independently.
#In diploid procumbens and webbiana, these correspond to blast searches that match only one
#contig from the pool assembly, as found in the xx* files.
#In tetraploid patellaris, there can be two or fewer contigs.
#The number of contigs present in any given blast result are counted by finding lines in the
#xx* files that start with 'jcf', the prefix given to contigs in the pool assemblies.

#tetraploids
for i in $(seq 50 1 52);
  do echo -n "$i ";
    s=0;
    mkdir "$i"/"$i"scogenes; #house a copy of the single copy orthologs here
    >"$i"/"$i"scos.txt;
    for j in $(ls "$i"/"$i"goodgenes/xx*);
      do a=$(grep ^jcf "$j" | wc -l);
        if (( a <= 2 ));
        then (( s = s + 1 ));
          echo "$j" >> "$i"/"$i"scos.txt;
          cp "$j" "$i"/"$i"scogenes;
        fi;
      done
    echo $s;
  done;

#diploids
for i in $(seq 53 1 55);
  do echo -n "$i ";
    s=0;
    mkdir "$i"/"$i"scogenes; #house a copy of the single copy orthologs here
    >"$i"/"$i"scos.txt;
    for j in $(ls "$i"/"$i"goodgenes/xx*);
      do a=$(grep ^jcf "$j" | wc -l);
        if (( a == 1 ));
        then (( s = s + 1 ));
          echo "$j" >> "$i"/"$i"scos.txt;
          cp "$j" "$i"/"$i"scogenes;
        fi;
      done
    echo $s;
  done;
			PoolID, Number of orthologous genes
			50 7289
			51 6031
			52 6721
			53 6953
			54 7797
			55 8617

			As expected there are many fewer single copy orthologs of EL10 genes in tetraploid P. patellaris (50,51,52)
			if you treat it as a diploid, since "single copy" genes are present in both genomes
			50 1518
			51 1002
			52 1254

#backup b100000yescp folder onto NAS at /share/Public/Data/PatReeves/PatellifoliaIlluminaAnalysis/reeves/FlashedReadArchive









### REMOVE BELOW BEFORE UPLOADING ###

###YOU MIGHT STILL NEED IT IF YOU GO ON TO QUANTIFY MAJOR ALLELE DIFFERENCES BETWEEN POOLS AT EACH GENE###

###URHERE use sgrep procedure below, you may not need this in the meta since having those reads
#might not have anything to do with the manuscript, but perform it anyway so you get the right scoreads.fa
#file in b100000yescp
#Use sgrep on sorted FASTA from FRA, *frasorted.fa to extract all reads matching orthologous genes
cd /scratch/reevesp;
rsync -aP admin@10.177.9.14:"/share/Public/Data/PatReeves/PatellifoliaIlluminaAnalysis/reeves/FlashedReadArchive/5*/5*fra/5*fraFinal/5*frasorted.fa" .;

cd /scratch/reevesp/patellifolia/blastdbflashedreads;
wd=$(pwd);
for i in $(seq 53 1 55);
  do echo "$i";
    > /scratch/reevesp/"$i"scoreads.fa;
    (grep ^rp "$wd"/"$i"/"$i"scogenes/xx* | cut -d' ' -f1 | cut -d: -f2 | sed 's/^/>/' \
         | parallel  'LC_ALL=C sgrep {} /scratch/reevesp/'"$i"'frasorted.fa | tr " " "\n"') \
         >> /scratch/reevesp/"$i"scoreads.fa;
  done;

#how many reads were extracted
cd /share/space/reevesp/patellifolia/blastdbflashedreads;
wd=$(pwd);
echo "pool readsexp readsfound";
for i in $(seq 50 1 55);
  do a=$(grep ^rp "$wd"/b100000yescp/"$i"/"$i"scogenes/xx* | wc -l); #number of reads expected from scogenes/xx*
    b=$(grep ^'>' "$wd"/b100000yescp/"$i"/"$i"scoreads.fa | wc -l); #number of reads found after extracting from *fra.fa
    echo "$i $a $b";
  done;
			pool readsexp readsfound
			50 261232 261232
			51 131744 131744
			52 199342 199342
			53 1821256 1821256
			54 1679330 1679330
			55 1842842 1842842

### REMOVE ABOVE BEFORE UPLOADING ###











### IDENTIFICATION OF SUGAR BEET EL10 GENOME HOMOLOGS IN PATELLIFOLIA ###

#calculate sequence length for nine EL10 v1.0 chromosomes
for i in $(seq 1 1 9);
  do echo -n Chr"$i ";
    grep -A1 ^'>'Chr$i Bvulgaris_548_EL10_1.0.fa | tail -1 | awk '{print length}';
  done;
			Chr1 58086001
			Chr2 54971872
			Chr3 54100447
			Chr4 61163185
			Chr5 59224585
			Chr6 65096967
			Chr7 57353724
			Chr8 57938902
			Chr9 52180088
			
#fragment chromosomes into consecutive 1000 bp fragments, takes about 4 minutes
time ( for i in $(seq 1 1 9);
  do echo Chr"$i";
    grep -A1 ^'>'Chr"$i" Bvulgaris_548_EL10_1.0.fa | tail -1 | sed -r 's/(.{1000})/\1\n/g' > Chr"$i".tmp;

    #calculate row labels
    a=$(wc -l Chr"$i".tmp | awk '{print $1}'); #number of lines for split chromosome
    b=$(grep ^'>'Chr"$i" Bvulgaris_548_EL10_1.0.fa | sed 's/>//'); #chromosome name
    c=1; #bp position of start of sequence fragment
    e="";
    for j in $(seq 1 1 $a);
      do d=">"$c"_$b"; #name of sequence fragment
        e+="$d"$'\n';#accumulate names for each sequence fragment
        c=$(( c + 1000 ));
      done;
    
    e=$(echo "$e" | sed 's/^$//'); #remove trailing newline
    paste -d' ' <(echo "$e") Chr"$i".tmp | tr " " "\n" > Chr"$i"_Bvulgaris_548_EL10_1.0.fa; #assemble output file
    rm Chr"$i".tmp;
  done;
)


#number of fragments is 1/1000x the number of bp shown above, so like 52-65K.
#Bvulgaris_548_EL10_1.0.cds_primaryTranscriptOnly.fa had 24K fragments of similar size so expect
#initial processing of genome to take 2x9chr = 18x times as long as transcripts (which took ~10 hrs for 6 pools and one query file)
#That means ~20 hrs if parallelized by pool within a node and by chromosome across nodes, or 180 hours if all Chromosomes run in series.
#Runs performed on ceres HPC using SLURM submit file as described below.


#prototyping for ceres
#use mem nodes with 80 cores, 666GB memory, all pools on one node, each Chr as a separate job, takes about 14 hours per job
#Use this slurm routine with blastdbs made from flashed reads 5*fra.fa.*

####BEGIN SLURM####
#!/bin/bash
#SBATCH --job-name="orthov1" #name of the job submitted
#SBATCH -p scavenger-mem,scavenger-mem768 #mem,mem768 #name of the queue you are submitting job to scavenger-mem,scavenger-mem768
#SBATCH -N 1 #number of nodes in this job
#SBATCH -n 80 #number of cores/tasks in this job
#SBATCH --mem 666G
#SBATCH -t 7-00:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=pat.reeves@ars.usda.gov #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "/dev/null" # standard out goes nowhere
#SBATCH -e "stderr.%j.%N.%A.%a" #optional but it prints out standard error

# submit this as an slurm array like: sbatch --array=1-9 orthov.slm
# submission like this induces a new variable SLURM_ARRAY_TASK_ID with values in the --array range.
# thereby, 9 slurm jobs are queued.  The variable can be used to specify folders with unique
# input files.  Each array job will occupy a different node.

#log to stderr
date >&2;

# stage data to /local/scratch ($TMPDIR)
mkdir $TMPDIR/blastdbflashedreads;
#copy blastdbs (about 8 minutes for 1 pool's blastdb, 42 minutes for all on n120 mem node, faster on new short nodes)
for i in $(seq 50 1 55);
  do rsync -aP /home/pat.reeves/patellifolia/OrthoVariant/blastdbflashedreads/"$i"frablastdb/"$i"fra.fa.* $TMPDIR/blastdbflashedreads/"$i"frablastdb;
  done;
#cp query files
rsync -aP /home/pat.reeves/patellifolia/OrthoVariant/blastdbflashedreads/Chr"$SLURM_ARRAY_TASK_ID"_Bvulgaris_548_EL10_1.0.fa $TMPDIR/blastdbflashedreads;

# load modules
module load parallel/20151022;
module load blast+/2.6.0;
 
# add regular job commands
g="Bvulgaris_548_EL10_1.0.fa"; #query file root
pqf="$TMPDIR/blastdbflashedreads"; #path to folder with query file

seq 50 1 55 | parallel 'blastn -num_threads 14 -num_descriptions 10000 -num_alignments 0 \
                                    -db '"$pqf"'/{}frablastdb/{}fra.fa \
                                    -query '"$pqf/Chr"$SLURM_ARRAY_TASK_ID"_$g"' \
                                    -out '"$pqf"'/{}fraorthov_'"Chr"$SLURM_ARRAY_TASK_ID"_$g"'.txt';

# copy output data off of local scratch and then remove all data from local scratch
rsync -aP "$pqf"/5[0-5]fraorthov_Chr"$SLURM_ARRAY_TASK_ID"_"$g".txt /home/pat.reeves/patellifolia/OrthoVariant;
  
#log to stderr
date >&2;

#End of file

####END SLURM####


#to "echo" the above slurm script to a file from the ceres terminal (instead of saving on local machine and rsyncing over)
cat > orthov.slm <<'EOF'
##PASTE SLURM SCRIPT###
EOF


#calculate total number of expected queries per chromosome
#on ceres:
cd /home/pat.reeves/patellifolia/OrthoVariant/blastdbflashedreads;
for i in $(seq 1 1 9);
  do echo -n "$i ";
    grep ^'>' Chr"$i"_Bvulgaris_548_EL10_1.0.fa | wc -l;
  done;
			Chromosome, NumQueries expected
			1 58087
			2 54972
			3 54101
			4 61164
			5 59225
			6 65097
			7 57354
			8 57939
			9 52181

#calculate number of queries actually conducted per chromosome, and number of queries that failed to return any hits
cd /home/pat.reeves/patellifolia/OrthoVariant;
(ls *fa.txt | parallel 'echo "{} "$(grep Query {} | wc -l)" "$(grep "No hits" {} | wc -l)') | sort -t_ -k2,2;

			All queries successful, note fairly substantial proportion (~25%) of queries
			are not alignable to any phased reads or pool assembly contigs present in the FRA.

			pool_chromosome, queries conducted, queries with no hits
			50fraorthov_Chr1_Bvulgaris_548_EL10_1.0.fa.txt 58087 14283
			51fraorthov_Chr1_Bvulgaris_548_EL10_1.0.fa.txt 58087 15935
			52fraorthov_Chr1_Bvulgaris_548_EL10_1.0.fa.txt 58087 14915
			53fraorthov_Chr1_Bvulgaris_548_EL10_1.0.fa.txt 58087 11305
			54fraorthov_Chr1_Bvulgaris_548_EL10_1.0.fa.txt 58087 15573
			55fraorthov_Chr1_Bvulgaris_548_EL10_1.0.fa.txt 58087 16055
			50fraorthov_Chr2_Bvulgaris_548_EL10_1.0.fa.txt 54972 13906
			51fraorthov_Chr2_Bvulgaris_548_EL10_1.0.fa.txt 54972 15332
			52fraorthov_Chr2_Bvulgaris_548_EL10_1.0.fa.txt 54972 14423
			53fraorthov_Chr2_Bvulgaris_548_EL10_1.0.fa.txt 54972 11187
			54fraorthov_Chr2_Bvulgaris_548_EL10_1.0.fa.txt 54972 14921
			55fraorthov_Chr2_Bvulgaris_548_EL10_1.0.fa.txt 54972 15451
			50fraorthov_Chr3_Bvulgaris_548_EL10_1.0.fa.txt 54101 14002
			51fraorthov_Chr3_Bvulgaris_548_EL10_1.0.fa.txt 54101 15400
			52fraorthov_Chr3_Bvulgaris_548_EL10_1.0.fa.txt 54101 14428
			53fraorthov_Chr3_Bvulgaris_548_EL10_1.0.fa.txt 54101 11301
			54fraorthov_Chr3_Bvulgaris_548_EL10_1.0.fa.txt 54101 14950
			55fraorthov_Chr3_Bvulgaris_548_EL10_1.0.fa.txt 54101 15380
			50fraorthov_Chr4_Bvulgaris_548_EL10_1.0.fa.txt 61164 15863
			51fraorthov_Chr4_Bvulgaris_548_EL10_1.0.fa.txt 61164 17677
			52fraorthov_Chr4_Bvulgaris_548_EL10_1.0.fa.txt 61164 16562
			53fraorthov_Chr4_Bvulgaris_548_EL10_1.0.fa.txt 61164 12674
			54fraorthov_Chr4_Bvulgaris_548_EL10_1.0.fa.txt 61164 17229
			55fraorthov_Chr4_Bvulgaris_548_EL10_1.0.fa.txt 61164 17708
			50fraorthov_Chr5_Bvulgaris_548_EL10_1.0.fa.txt 59225 15042
			51fraorthov_Chr5_Bvulgaris_548_EL10_1.0.fa.txt 59225 16739
			52fraorthov_Chr5_Bvulgaris_548_EL10_1.0.fa.txt 59225 15665
			53fraorthov_Chr5_Bvulgaris_548_EL10_1.0.fa.txt 59225 12022
			54fraorthov_Chr5_Bvulgaris_548_EL10_1.0.fa.txt 59225 16361
			55fraorthov_Chr5_Bvulgaris_548_EL10_1.0.fa.txt 59225 16932
			50fraorthov_Chr6_Bvulgaris_548_EL10_1.0.fa.txt 65097 16583
			51fraorthov_Chr6_Bvulgaris_548_EL10_1.0.fa.txt 65097 18478
			52fraorthov_Chr6_Bvulgaris_548_EL10_1.0.fa.txt 65097 17270
			53fraorthov_Chr6_Bvulgaris_548_EL10_1.0.fa.txt 65097 13127
			54fraorthov_Chr6_Bvulgaris_548_EL10_1.0.fa.txt 65097 18023
			55fraorthov_Chr6_Bvulgaris_548_EL10_1.0.fa.txt 65097 18513
			50fraorthov_Chr7_Bvulgaris_548_EL10_1.0.fa.txt 57354 14590
			51fraorthov_Chr7_Bvulgaris_548_EL10_1.0.fa.txt 57354 16140
			52fraorthov_Chr7_Bvulgaris_548_EL10_1.0.fa.txt 57354 15134
			53fraorthov_Chr7_Bvulgaris_548_EL10_1.0.fa.txt 57354 11560
			54fraorthov_Chr7_Bvulgaris_548_EL10_1.0.fa.txt 57354 15665
			55fraorthov_Chr7_Bvulgaris_548_EL10_1.0.fa.txt 57354 16139
			50fraorthov_Chr8_Bvulgaris_548_EL10_1.0.fa.txt 57939 14430
			51fraorthov_Chr8_Bvulgaris_548_EL10_1.0.fa.txt 57939 16079
			52fraorthov_Chr8_Bvulgaris_548_EL10_1.0.fa.txt 57939 14967
			53fraorthov_Chr8_Bvulgaris_548_EL10_1.0.fa.txt 57939 11536
			54fraorthov_Chr8_Bvulgaris_548_EL10_1.0.fa.txt 57939 15656
			55fraorthov_Chr8_Bvulgaris_548_EL10_1.0.fa.txt 57939 16178
			50fraorthov_Chr9_Bvulgaris_548_EL10_1.0.fa.txt 52181 13622
			51fraorthov_Chr9_Bvulgaris_548_EL10_1.0.fa.txt 52181 14984
			52fraorthov_Chr9_Bvulgaris_548_EL10_1.0.fa.txt 52181 13995
			53fraorthov_Chr9_Bvulgaris_548_EL10_1.0.fa.txt 52181 11069
			54fraorthov_Chr9_Bvulgaris_548_EL10_1.0.fa.txt 52181 14737
			55fraorthov_Chr9_Bvulgaris_548_EL10_1.0.fa.txt 52181 15115





### FILTER BLASTN RESULTS TO FIND SINGLE COPY ORTHOLOG GENOMIC REGIONS BETWEEN EL10 AND PATELLIFOLIA ###

#find EL10 1kb fragments for which there is only 1 patellifolia contig that matches (2 or
#fewer in tetraploids). These are defined as "single copy" orthologs
#takes a couple minutes 

#tetraploids
mypp() {
        i=$1; #input file name
        p=$(echo "$i" | cut -c1-2); #pool name
        c=$(echo "$i" | cut -d_ -f2); #chr name
        awk 'BEGIN {RS="Query= "} {print $1,gsub(/\njcf/,"")}' "$i" | grep " "[12]$ > "$p"_"$c"_scos.txt;
        echo "$i done";
}
export -f mypp;

cd /home/pat.reeves/patellifolia/OrthoVariant/blastnresultsflashedreads;
screen; sscavenger; #get 72 cores on ceres
ls 5[0-2]fraorthov_Chr[0-9]_Bvulgaris_548_EL10_1.0.fa.txt | parallel --env mypp mypp;

#diploids
mypp() {
        i=$1; #input file name
        p=$(echo "$i" | cut -c1-2); #pool name
        c=$(echo "$i" | cut -d_ -f2); #chr name
        awk 'BEGIN {RS="Query= "} {print $1,gsub(/\njcf/,"")}' "$i" | grep " 1$" > "$p"_"$c"_scos.txt;
        echo "$i done";
}
export -f mypp;

cd /home/pat.reeves/patellifolia/OrthoVariant/blastnresultsflashedreads;
screen; sscavenger; #get 72 cores on ceres
ls 5[3-5]fraorthov_Chr[0-9]_Bvulgaris_548_EL10_1.0.fa.txt | parallel --env mypp mypp;

#count "single copy" orthologous regions
for i in *scos.txt; do wc -l "$i"; done; #how many scos per pool per chr?
			NumOrthologousRegions, Pool_Chromosome
			12248 50_Chr1_scos.txt
			11749 50_Chr2_scos.txt
			11998 50_Chr3_scos.txt
			13344 50_Chr4_scos.txt
			12884 50_Chr5_scos.txt
			13824 50_Chr6_scos.txt
			12135 50_Chr7_scos.txt
			12285 50_Chr8_scos.txt
			11192 50_Chr9_scos.txt
			11972 51_Chr1_scos.txt
			11638 51_Chr2_scos.txt
			11822 51_Chr3_scos.txt
			12946 51_Chr4_scos.txt
			12690 51_Chr5_scos.txt
			13607 51_Chr6_scos.txt
			11904 51_Chr7_scos.txt
			12125 51_Chr8_scos.txt
			11019 51_Chr9_scos.txt
			12115 52_Chr1_scos.txt
			11693 52_Chr2_scos.txt
			12031 52_Chr3_scos.txt
			13335 52_Chr4_scos.txt
			12666 52_Chr5_scos.txt
			13708 52_Chr6_scos.txt
			12029 52_Chr7_scos.txt
			12450 52_Chr8_scos.txt
			11158 52_Chr9_scos.txt
			8586 53_Chr1_scos.txt
			8330 53_Chr2_scos.txt
			8302 53_Chr3_scos.txt
			9414 53_Chr4_scos.txt
			8901 53_Chr5_scos.txt
			9636 53_Chr6_scos.txt
			8459 53_Chr7_scos.txt
			8545 53_Chr8_scos.txt
			7969 53_Chr9_scos.txt
			8983 54_Chr1_scos.txt
			8621 54_Chr2_scos.txt
			8788 54_Chr3_scos.txt
			9827 54_Chr4_scos.txt
			9423 54_Chr5_scos.txt
			10044 54_Chr6_scos.txt
			8896 54_Chr7_scos.txt
			9029 54_Chr8_scos.txt
			8281 54_Chr9_scos.txt
			9368 55_Chr1_scos.txt
			8959 55_Chr2_scos.txt
			9233 55_Chr3_scos.txt
			10257 55_Chr4_scos.txt
			9940 55_Chr5_scos.txt
			10595 55_Chr6_scos.txt
			9296 55_Chr7_scos.txt
			9515 55_Chr8_scos.txt
			8786 55_Chr9_scos.txt

			fewer for tetraploids if treated as diploids, as expected (5-7K vs 8-10K)
			6217 50_Chr1_scos.txt
			5822 50_Chr2_scos.txt
			6038 50_Chr3_scos.txt
			6892 50_Chr4_scos.txt
			6500 50_Chr5_scos.txt
			6933 50_Chr6_scos.txt
			6208 50_Chr7_scos.txt
			6195 50_Chr8_scos.txt
			5438 50_Chr9_scos.txt
			6126 51_Chr1_scos.txt
			5899 51_Chr2_scos.txt
			6054 51_Chr3_scos.txt
			6705 51_Chr4_scos.txt
			6383 51_Chr5_scos.txt
			6953 51_Chr6_scos.txt
			6251 51_Chr7_scos.txt
			6145 51_Chr8_scos.txt
			5498 51_Chr9_scos.txt
			6225 52_Chr1_scos.txt
			5882 52_Chr2_scos.txt
			6205 52_Chr3_scos.txt
			6846 52_Chr4_scos.txt
			6564 52_Chr5_scos.txt
			7004 52_Chr6_scos.txt
			6224 52_Chr7_scos.txt
			6493 52_Chr8_scos.txt
			5640 52_Chr9_scos.txt

#Make a plottable file showing regions of EL10 with scos in Patellifolia, ~4 min
#Get a "depth" (phased read count) for each pool, at each 1kb position in EL10 genome. ~3 min

mypp() {
       j=$1; #input file like 5[0-5]fraorthov_Chr[0-9]_Bvulgaris_548_EL10_1.0.fa.txt
       i=$(echo "$j" | cut -c1-2); #pool number
       c=$(echo "$j" | cut -d_ -f2 | sed 's/Chr//'); #chromosome number
       
       awk 'BEGIN {RS="Query= "} {print $1 " " gsub(/\nrp/,"")}' "$j" \
         | tail -n +2 | sed 's/_Chr'"$c"'_EL10_PGA_scaffold[0-9] /:/g' > "$i".Chr"$c".2.tmp;
}
export -f mypp

cd /home/pat.reeves/patellifolia/OrthoVariant/blastnresultsflashedreads;
time ls 5[0-5]fraorthov_Chr[0-9]_Bvulgaris_548_EL10_1.0.fa.txt | parallel --env mypp mypp;
#assemble awk output into a plottable file for each chromosome
for i in {1..9};
  do paste -d$'\t' <(awk -F: '{print $1}' 50.Chr"$i".2.tmp) \
                   <(awk -F: '{print $2}' 50.Chr"$i".2.tmp) <(awk -F: '{print $2}' 51.Chr"$i".2.tmp) <(awk -F: '{print $2}' 52.Chr"$i".2.tmp) \
                   <(awk -F: '{print $2}' 53.Chr"$i".2.tmp) <(awk -F: '{print $2}' 54.Chr"$i".2.tmp) <(awk -F: '{print $2}' 55.Chr"$i".2.tmp) \
                   > Chr"$i"plot2.txt;
  done;

rm 5[0-5].Chr[0-9].2.tmp; #clean up

#reduce flashed read counts per region to scos only ~2min
mypp() {
       i=$1; #chromosome number
       a=$(cut -d$'\t' -f1 Chr"$i"plot.txt); #get the sco regions from the first plot files
       (echo "$a" | parallel --keep-order 'LC_ALL=C grep ^{}'"$'\t'"' Chr'"$i"'plot2.txt') > Chr"$i"plot3.txt;
}
export -f mypp;

cd /home/pat.reeves/patellifolia/OrthoVariant/blastnresultsflashedreads;
time seq 1 1 9 | parallel --env mypp mypp;


#back up work to nas


#on local machine, use R to make an animated bar plot sliding along the chromosome,
#showing phased read depth on y axis
cd /Users/wichita/Desktop/telework/patellifolia/OrthoVariantMap/plots/animations;


### BEGIN R ### #Log scale Y axis

options(error = recover)
rm(list=ls()) 

library(doParallel)
library(foreach)
setwd("/Users/wichita/Desktop/telework/patellifolia/OrthoVariantMap/plots")

for (c in seq(1,9,1))
{
  infilename=paste("Chr",c,"plot3.txt",sep="") #plot common scos only
#  infilename=paste("Chr",c,"plot2.txt",sep="") #plot depth at all 1kb fragments
  
  g=read.table(infilename, header=FALSE, sep="\t")
  colnames(g) = c("pos","a50","a51","a52","a53","a54","a55")
  
  #set up for parallel
  cores=detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterExport(cl,c("g"),envir=environment());
  clusterEvalQ(cl,{
    library(ggplot2)
    library(grid)
    library(gridExtra)
    })

#  for (i in seq(1,max(g$pos),10000)) #this is the step size in bp between each slide.
  foreach(i=seq(1,max(g$pos),10000)) %dopar%
  {
    gsub=g[g$pos>=i & g$pos<i+1e6,] #defines the number of 1kb regions to include in each slide, here include 1e6 bp worth of 1kb regions
    outfilename=paste(formatC(i,width=8,format="d",flag="0"),".jpg",sep="") #produce zero-padded outfile for proper sorting with ffmpeg
  
    xbreaks=round(seq(min(gsub$pos), max(gsub$pos), by = 500000)-1)
    
    #plot for pool 50
    ggp50=ggplot(data=gsub,aes(x=pos,y=a50)) + 
      geom_col()+
      ylab("patellaris") +
      xlab(paste("Chr",c,"position")) +
      scale_y_log10() +
      coord_cartesian(ylim=c(1,10000)) +
      scale_x_continuous(breaks=xbreaks)
  
    #plot for pool 51
    ggp51=ggplot(data=gsub,aes(x=pos,y=a51)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("patellaris") +
      scale_y_log10() +
      coord_cartesian(ylim=c(1,10000)) +
      scale_x_continuous(breaks=xbreaks)
  
    #plot for pool 52
    ggp52=ggplot(data=gsub,aes(x=pos,y=a52)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("patellaris") +
      scale_y_log10() +
      coord_cartesian(ylim=c(1,10000)) +
      scale_x_continuous(breaks=xbreaks)
  
    #plot for pool 53
    ggp53=ggplot(data=gsub,aes(x=pos,y=a53)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("procumbens") +
      scale_y_log10() +
      coord_cartesian(ylim=c(1,10000)) +
      scale_x_continuous(breaks=xbreaks)
  
    #plot for pool 54
    ggp54=ggplot(data=gsub,aes(x=pos,y=a54)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("procumbens") +
      scale_y_log10() +
      coord_cartesian(ylim=c(1,10000)) +
      scale_x_continuous(breaks=xbreaks)
  
    #plot for pool 55
    ggp55=ggplot(data=gsub,aes(x=pos,y=a55)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("webbiana") +
      scale_y_log10() +
      coord_cartesian(ylim=c(1,10000)) +
      scale_x_continuous(breaks=xbreaks)
  
    #arrange plots and print to jpeg
    jpeg(filename=paste(getwd(),"/animations/",outfilename, sep=""), width=7, height=7, units="in", res=144, bg="white")
    grid.draw(rbind(ggplotGrob(ggp55),ggplotGrob(ggp54),ggplotGrob(ggp53),ggplotGrob(ggp52),ggplotGrob(ggp51),ggplotGrob(ggp50),size="last"))
 
    dev.off()
    
  } #make slides
  
  #issue system calls to create animation, then clean up
 x=paste("ffmpeg -r 60 -pattern_type glob -i 'animations/*.jpg' -vcodec libx264 -pix_fmt yuv420p animations/Chr",c,"scodepthlog.mp4",sep="") #plot scos only
#  x=paste("ffmpeg -r 60 -pattern_type glob -i 'animations/*.jpg' -vcodec libx264 -pix_fmt yuv420p animations/Chr",c,"scodepth.mp4",sep="") #plot scos only
#  x=paste("ffmpeg -r 60 -pattern_type glob -i 'animations/*.jpg' -vcodec libx264 -pix_fmt yuv420p animations/Chr",c,"alldepth.mp4",sep="") #plot all 1kb fragments
  system(x)
  system("rm animations/*.jpg")
  
  stopCluster(cl)

} #iterate across chromosomes

### END R ### #Log scale Y axis


### BEGIN R ###  #Normal scale Y axis
options(error = recover)
rm(list=ls()) 

library(doParallel)
library(foreach)
setwd("/Users/wichita/Desktop/telework/patellifolia/OrthoVariantMap/plots")

for (c in seq(1,9,1))
{
  infilename=paste("Chr",c,"plot3.txt",sep="") #plot common scos only
#  infilename=paste("Chr",c,"plot2.txt",sep="") #plot depth at all 1kb fragments
  
  g=read.table(infilename, header=FALSE, sep="\t")
  colnames(g) = c("pos","a50","a51","a52","a53","a54","a55")
  
  #set up for parallel
  cores=detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterExport(cl,c("g"),envir=environment());
  clusterEvalQ(cl,{
    library(ggplot2)
    library(grid)
    library(gridExtra)
    })

#  for (i in seq(1,max(g$pos),10000)) #this is the step size in bp between each slide.
  foreach(i=seq(1,max(g$pos),10000)) %dopar%
  {
    gsub=g[g$pos>=i & g$pos<i+1e6,] #defines the number of 1kb regions to include in each slide, here include 1e6 bp worth of 1kb regions
    outfilename=paste(formatC(i,width=8,format="d",flag="0"),".jpg",sep="") #produce zero-padded outfile for proper sorting with ffmpeg
  
    #plot for pool 50
    ggp50=ggplot(data=gsub,aes(x=pos,y=a50)) + 
      geom_col()+
      ylab("patellaris") +
      xlab(paste("Chr",c,"position")) +
      coord_cartesian(ylim=c(0,1000)) +
      scale_x_continuous(breaks=round(seq(min(gsub$pos), max(gsub$pos), by = 500000)-1) )
  
    #plot for pool 51
    ggp51=ggplot(data=gsub,aes(x=pos,y=a51)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("patellaris") +
      coord_cartesian(ylim=c(0,1000)) +
      scale_x_continuous(breaks=round(seq(min(gsub$pos), max(gsub$pos), by = 500000)-1) )
  
    #plot for pool 52
    ggp52=ggplot(data=gsub,aes(x=pos,y=a52)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("patellaris") +
      coord_cartesian(ylim=c(0,1000)) +
      scale_x_continuous(breaks=round(seq(min(gsub$pos), max(gsub$pos), by = 500000)-1) )
  
    #plot for pool 53
    ggp53=ggplot(data=gsub,aes(x=pos,y=a53)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("procumbens") +
      coord_cartesian(ylim=c(0,1000)) +
      scale_x_continuous(breaks=round(seq(min(gsub$pos), max(gsub$pos), by = 500000)-1) )
  
    #plot for pool 54
    ggp54=ggplot(data=gsub,aes(x=pos,y=a54)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("procumbens") +
      coord_cartesian(ylim=c(0,1000)) +
      scale_x_continuous(breaks=round(seq(min(gsub$pos), max(gsub$pos), by = 500000)-1) )
  
    #plot for pool 55
    ggp55=ggplot(data=gsub,aes(x=pos,y=a55)) + 
      geom_col()+
      theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +
      ylab("webbiana") +
      coord_cartesian(ylim=c(0,1000)) +
      scale_x_continuous(breaks=round(seq(min(gsub$pos), max(gsub$pos), by = 500000)-1) )
  
    #arrange plots and print to jpeg
    jpeg(filename=paste(getwd(),"/animations/",outfilename, sep=""), width=7, height=7, units="in", res=144, bg="white")
    grid.draw(rbind(ggplotGrob(ggp55),ggplotGrob(ggp54),ggplotGrob(ggp53),ggplotGrob(ggp52),ggplotGrob(ggp51),ggplotGrob(ggp50),size="last"))
 
    dev.off()
    
  } #make slides
  
  #issue system calls to create animation, then clean up
  x=paste("ffmpeg -r 60 -pattern_type glob -i 'animations/*.jpg' -vcodec libx264 -pix_fmt yuv420p animations/Chr",c,"scodepth.mp4",sep="") #plot scos only
#  x=paste("ffmpeg -r 60 -pattern_type glob -i 'animations/*.jpg' -vcodec libx264 -pix_fmt yuv420p animations/Chr",c,"alldepth.mp4",sep="") #plot all 1kb fragments
  system(x)
  system("rm animations/*.jpg")
  
  stopCluster(cl)

} #iterate across chromosomes

### END R ### #Normal scale Y axis


#create a mosaic of all chromosome depth movies, log scale Y axis, all homologous regions
cd /Users/wichita/Desktop/telework/patellifolia/OrthoVariantMap/plots/animations 
ffmpeg \
   -i Chr1alldepth.mp4 \
   -i Chr2alldepth.mp4 \
   -i Chr3alldepth.mp4 \
   -i Chr4alldepth.mp4 \
   -i Chr5alldepth.mp4 \
   -i Chr6alldepth.mp4 \
   -i Chr7alldepth.mp4 \
   -i Chr8alldepth.mp4 \
   -i Chr9alldepth.mp4 \
   -filter_complex "[0:v][1:v][2:v][3:v][4:v][5:v][6:v][7:v][8:v]xstack=inputs=9:layout=0_0|w0_0|w0+w1_0|0_h0|w0_h0|w0+w1_h0|0_h0+h1|w0_h0+h1|w0+w1_h0+h1[v]" \
       -map "[v]" alldepthMosaic.mp4;

#create a mosaic of all chromosome depth movies, normal scale Y axis, single copy regions only
cd /Users/wichita/Desktop/telework/patellifolia/OrthoVariantMap/plots/animations 
ffmpeg \
   -i Chr1scodepth.mp4 \
   -i Chr2scodepth.mp4 \
   -i Chr3scodepth.mp4 \
   -i Chr4scodepth.mp4 \
   -i Chr5scodepth.mp4 \
   -i Chr6scodepth.mp4 \
   -i Chr7scodepth.mp4 \
   -i Chr8scodepth.mp4 \
   -i Chr9scodepth.mp4 \
   -filter_complex "[0:v][1:v][2:v][3:v][4:v][5:v][6:v][7:v][8:v]xstack=inputs=9:layout=0_0|w0_0|w0+w1_0|0_h0|w0_h0|w0+w1_h0|0_h0+h1|w0_h0+h1|w0+w1_h0+h1[v]" \
       -map "[v]" scodepthMosaic.mp4;

#create a mosaic of all chromosome depth movies, log scale Y axis, single copy regions only
cd /Users/wichita/Desktop/telework/patellifolia/OrthoVariantMap/plots/animations 
ffmpeg \
   -i Chr1scodepthlog.mp4 \
   -i Chr2scodepthlog.mp4 \
   -i Chr3scodepthlog.mp4 \
   -i Chr4scodepthlog.mp4 \
   -i Chr5scodepthlog.mp4 \
   -i Chr6scodepthlog.mp4 \
   -i Chr7scodepthlog.mp4 \
   -i Chr8scodepthlog.mp4 \
   -i Chr9scodepthlog.mp4 \
   -filter_complex "[0:v][1:v][2:v][3:v][4:v][5:v][6:v][7:v][8:v]xstack=inputs=9:layout=0_0|w0_0|w0+w1_0|0_h0|w0_h0|w0+w1_h0|0_h0+h1|w0_h0+h1|w0+w1_h0+h1[v]" \
       -map "[v]" scodepthlogMosaic.mp4;



### DETERMINE LOCATIONS OF HS1PRO1 AND HS4 AND BvBTC1 IN EL10 GENOME###

module load blast+/2.9.0;

#Hs1pro1
blastn -db /home/pat.reeves/patellifolia/EL10BlastDBs/1kb_Bvulgaris_548_EL10_1.0.fa -query /home/pat.reeves/patellifolia/seq/Hs1pro-1.fa -out hs1pro1.EL10.out.txt
			Query= U79733.1 Beta procumbens nematode resistance (Hs1pro-1) mRNA,
			complete cds

			Length=1450
																				  Score        E
			Sequences producing significant alignments:                          (Bits)     Value

			53555001_Chr2_EL10_PGA_scaffold6                                      523        2e-146
			53554001_Chr2_EL10_PGA_scaffold6                                      243        5e-62 


#putative Hs4
blastn -db /home/pat.reeves/patellifolia/EL10BlastDBs/1kb_Bvulgaris_548_EL10_1.0.fa -query /home/pat.reeves/patellifolia/seq/XM_010669575.fa -out hs4putative.EL10.out.txt
			Query= XM_010669575.1 PREDICTED: Beta vulgaris subsp. vulgaris hypothetical
			protein (LOC104884871), mRNA

			Length=1292
																				  Score        E
			Sequences producing significant alignments:                          (Bits)     Value

			44495001_Chr2_EL10_PGA_scaffold6                                      935        0.0   
			44491001_Chr2_EL10_PGA_scaffold6                                      418        7e-115
			44494001_Chr2_EL10_PGA_scaffold6                                      311        1e-82 
			44493001_Chr2_EL10_PGA_scaffold6                                      237        2e-60 
			44490001_Chr2_EL10_PGA_scaffold6                                      117        3e-24 

#BvBTC1 Bolting gene: 32609001_Chr2_EL10_PGA_scaffold6:32619001_Chr2_EL10_PGA_scaffold6
blastn -db /home/pat.reeves/patellifolia/EL10BlastDBs/1kb_Bvulgaris_548_EL10_1.0.fa -query /home/pat.reeves/patellifolia/seq/BvBTC1.fa -out BvBTC1.EL10.out.txt
			Query= HQ709091.1 Beta vulgaris subsp. vulgaris genotype KWS2320 bolting
			time control 1 (BTC1) gene, complete cds

			Length=14333
																				  Score        E
			Sequences producing significant alignments:                          (Bits)     Value

			32615001_Chr2_EL10_PGA_scaffold6                                      1834       0.0   
			32616001_Chr2_EL10_PGA_scaffold6                                      1825       0.0   
			32613001_Chr2_EL10_PGA_scaffold6                                      1821       0.0   
			32611001_Chr2_EL10_PGA_scaffold6                                      1820       0.0   
			32612001_Chr2_EL10_PGA_scaffold6                                      1808       0.0   
			32618001_Chr2_EL10_PGA_scaffold6                                      1788       0.0   
			32609001_Chr2_EL10_PGA_scaffold6                                      1768       0.0   
			32610001_Chr2_EL10_PGA_scaffold6                                      1755       0.0   
			32617001_Chr2_EL10_PGA_scaffold6                                      1659       0.0   
			32619001_Chr2_EL10_PGA_scaffold6                                      1199       0.0   
blastn -db /home/pat.reeves/patellifolia/EL10BlastDBs/1kb_Bvulgaris_548_EL10_1.0.fa -query /home/pat.reeves/patellifolia/seq/BvBTC1mRNA.fa -out BvBTC1mRNA.EL10.out.txt
			Query= HQ709091.1:2001-2131,2259-2297,2693-3192,4100-4261,4421-4554,
			4851-5012,10105-10285,10367-10802,10877-11699,11991-12338 Beta
			vulgaris subsp. vulgaris genotype KWS2320 bolting time control 1
			(BTC1) gene, complete cds

			Length=2916
																				  Score        E
			Sequences producing significant alignments:                          (Bits)     Value

			32611001_Chr2_EL10_PGA_scaffold6                                      1445       0.0   
			32618001_Chr2_EL10_PGA_scaffold6                                      915        0.0   
			32612001_Chr2_EL10_PGA_scaffold6                                      802        0.0   
			32610001_Chr2_EL10_PGA_scaffold6                                      625        9e-177
			32617001_Chr2_EL10_PGA_scaffold6                                      303        5e-80 
			32616001_Chr2_EL10_PGA_scaffold6                                      294        3e-77 
			32619001_Chr2_EL10_PGA_scaffold6                                      231        2e-58 
			49598001_Chr1_EL10_PGA_scaffold3                                      231        2e-58 
			49599001_Chr1_EL10_PGA_scaffold3                                      211        3e-52 



### HS1PRO1 ALLELE MINING ###
### DETERMINE CONTIGS CONTAINING HS1PRO1 IN PATELLIFOLIA POOL ASSEMBLIES ###

module load blast+/2.9.0;

cd /home/pat.reeves/patellifolia/FlashedReadArchive;

#Hs1pro1
for i in $(seq 50 1 55);
  do blastn -db /home/pat.reeves/patellifolia/FlashedReadArchive/"$i"fraFinal/"$i"Hs1pro1REVdref.fasta -query /home/pat.reeves/patellifolia/seq/Hs1pro-1.fa -out "$i"mapxHs1pro1.txt;
  done;

#consolidated output from blast
#all searches above return 1 hit in procumbens/webbiana, two hits in all three patellaris because it is tetraploid
			Query= U79733.1 Beta procumbens nematode resistance (Hs1pro-1) mRNA,
			complete cds
			Length=1450
																				  Score        E
			Sequences producing significant alignments:                          (Bits)     Value
			50jcf7180007932120rev                                                 2451       0.0  
			50jcf7180007994428rev                                                 2039       0.0  
			51jcf7180007742276                                                    2451       0.0  
			51jcf7180007654607                                                    2045       0.0  
			52jcf7180008015606rev                                                 2423       0.0  
			52jcf7180007917886                                                    2045       0.0  
			53jcf7180008859096rev                                                 2495       0.0  
			54jcf7180009178902rev                                                 2446       0.0  
			55jcf7180008573862rev                                                 2446       0.0  


#Manual inspection of aligned sequences using sequencher showed shared indels can be used to
#distinguish orthologs.
#Within tetraploid patellaris one of the two sequences returned by blast is orthologous to
#the single sequence found in procumbens/webbiana

			Hs1pro1 locus 1 orthologs, and spans containing the Hs1pro1L1 gene are:
			50jcf7180007932120rev_1493-5327
			51jcf7180007742276_10788-14621
			52jcf7180008015606rev_7103-10936
			53jcf7180008859096rev_1407-5320
			54jcf7180009178902rev_15333-19226
			55jcf7180008573862rev_2380-6224

			Hs1pro1 locus 2 orthologs, and spans containing Hs1pro1L2 gene (not used here) are:
			50jcf7180007994428rev_6873-10755
			51jcf7180007654607_6772-10654
			52jcf7180007917886_1297-5179
			
#After identifying which contigs contain Hs1pro1 locus 1 in each pool assembly, and the
#coordinates within those contigs that contain the Hs1pro1L1 gene, extract the phased reads
#that map to that region from the FRA map. Just use samtools for this since the reads are
#already phased and quality filtered by virtue of being in the map		
			
cd /home/pat.reeves/patellifolia/AlleleMining/Hs1pro1;
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/50fraFinal/50fra.bam "50jcf7180007932120rev:1493-5327" > 50Hs1pro1L1.bam 
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/51fraFinal/51fra.bam "51jcf7180007742276:10788-14621" > 51Hs1pro1L1.bam 
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/52fraFinal/52fra.bam "52jcf7180008015606rev:7103-10936" > 52Hs1pro1L1.bam 
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/53fraFinal/53fra.bam "53jcf7180008859096rev:1407-5320" > 53Hs1pro1L1.bam 
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/54fraFinal/54fra.bam "54jcf7180009178902rev:15333-19226" > 54Hs1pro1L1.bam 
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/55fraFinal/55fra.bam "55jcf7180008573862rev:2380-6224" > 55Hs1pro1L1.bam 

#extract phased reads from bam files into a single fasta file
seq 50 1 55 | parallel --keep-order 'echo {}; samtools fasta {}Hs1pro1L1.bam > {}Hs1pro1L1.fa';
			50
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 822 reads
			51
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 425 reads
			52
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 645 reads
			53
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 1345 reads
			54
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 914 reads
			55
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 1119 reads
			
#map these phased reads to a single reference gene sequence (choose the longest one)
#for Hs1pro1L1 use 51jcf7180007742276_ref.txt, indexed with bwa index
cd /home/pat.reeves/patellifolia/AlleleMining/Hs1pro1;
mkdir bam;
pd=$(pwd);
t=40; #threads
for j in 50Hs1pro1L1.fa 51Hs1pro1L1.fa 52Hs1pro1L1.fa 53Hs1pro1L1.fa 54Hs1pro1L1.fa 55Hs1pro1L1.fa;
    do rg=$(cut -c1-2 <<<"$j"); #get the read group name 5[0-5]
      echo "$rg";
      readgroup="@RG\tID:$rg\tSM:$rg\tPL:illumina\tLB:na\tPU:na";
      bwa mem -R "$readgroup" -t "$t" "$pd"/ref/51jcf7180007742276_ref.txt "$pd"/"$j" | \
                  samtools sort -O BAM --threads "$t" | \
                  samtools view -F 2048 -O BAM > "$pd"/bam/"$rg"Hs1pro1l1.finalaln.bam.TMP; #-F 2048 excludes supplementary alignments
                  samtools index "$pd"/bam/"$rg"Hs1pro1l1.finalaln.bam.TMP;
    done;
sambamba merge -t40 -p "$pd"/bam/Hs1pro1l1.finalaln.merge.bam "$pd"/bam/5[0-5]Hs1pro1l1.finalaln.bam.TMP;

### END DETERMINE CONTIGS CONTAINING HS1PRO1 IN PATELLIFOLIA POOL ASSEMBLIES ###



### OPTIONAL: COUNT TOTAL, DISTINCT, AND PRIVATE ALLELES ###
#rsync results from HPC ceres to HPC blip so you can get access to multi-node gnu parallel

#run a final hapx to count microhaplotypes at each position (takes ~35 minutes)
cd /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1; #go here just to get the working directory specified
pd=$(pwd);
time hapx.sh -r "$pd"/ref/51jcf7180007742276_ref.txt \
        -b "$pd"/bam/Hs1pro1l1.finalaln.merge.bam \
        -o Hs1pro1l1Finalhapx \
        -f 0 -q 1 -x -ssh /home/reevesp/machines \
        -s <(for i in $(seq 9000 1 16000); do echo 51jcf7180007742276:"$i"-"$i"; done;);


#postprocess haploblock counts in log.txt into something plottable
cd /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1/Hs1pro1l1Finalhapx;

### BEGIN BASH ###
mytd1() {
        i=$1; #a position on the reference is incoming
        a=$(grep "$i"\\."$k" "$pd"/numblocks.txt | cut -d: -f4); #number of unique haploblock read pairs mapped to position i in the readgroup
        if [[ "$a" == "" ]]; then a="?"; fi; #if no data at position i set number of haploblock read pairs to ?
        echo "$i"."$k $a"; #report result to parallel statement
}
export -f mytd1;

mytd2() {
        i=$1; #a position on the reference is incoming
        b=$(grep "$i"\\."$k" "$pd"/numblocks.txt | cut -d$'\t' -f2 | cut -d: -f1); #total number haploblock read pairs mapped to position i in the readgroup
        if [[ "$b" == "" ]]; then b="?"; fi; #if no data at position i set to ?
        echo "$i"."$k $b"; #report result to parallel statement
}
export -f mytd2;

mypa() {
       i=$1; #contig:site-range
       j=$(grep "$i" "$pd"/counts.txt | grep -v global);
       names=$(cut -d$'\t' -f1 <<<"$j" | sed 's/$/ 0/');
       counts=$(cut -d$'\t' -f2 <<<"$j");
       nc=$(head -1 <<<"$counts" | awk -F: '{print NF}'); #number of alleles
       nr=$(wc -l <<<"$counts"); #number of read groups
       nz=$(( $nr - 1 )); #number of zeroes needed in a column of allele counts for a private allele to exist
       for m in $(seq 1 1 $nc);
         do col=$(cut -d: -f$m <<<"$counts"); #extract the column of allele counts
           fl=$(grep -n -v 0 <<<"$col"); #show line numbers where alleles are present
           
           #if only one line of the column of data has alleles, it is a private allele
           if [[ $(echo "$fl" | wc -l) == 1 ]];
           then r=$(cut -d: -f1 <<<"$fl"); #row number of private allele
             names=$(awk -F' ' -v r=$r '{if (NR==r) $2++; print}' <<<"$names"); #index up by one the row with the private allele
           fi;
         done;
         
       #report result to parallel statement
       echo "$names";
}
export -f mypa;

mypadpa() {
          aa=$1;
          cc=$(grep "$aa.$dd" "$pd"/names.txt);
          if [[ "$cc" == "" ]];
          then echo "$aa.$dd ?"; #if readgroup not found print a ?
          else echo "$cc";
          fi;
}
export -f mypadpa;


pd=$(pwd); export pd;

#count total and distinct haploblocks
grep ^# log.txt | tr '-' '_' > numblocks.txt;
a=$(rev numblocks.txt | cut -d. -f2- | rev | uniq); #list of contig:site-ranges
b=$(rev numblocks.txt | cut -d. -f1 | rev | cut -d$'\t' -f1 | sort -u); export b; #determine set of possible read groups


#iterate over read groups to count haploblocks
for j in $b;
  do k="$j"; export k; #do this so iterator can be sent to nodes using --sshloginfile
    echo "$j: counting distinct alleles";
    >"$pd"/"$j".distallel.txt; #initialize output file with a header for count of distinct alleles within populations
    echo "$a" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env mytd1 --env pd --env k /home/reevesp/bin/parallel --jobs 96 --env mytd1 --env pd --env k mytd1 >> "$pd"/"$j".distallel.txt;
    #echo "$a" | parallel --jobs 1 --pipe -N960 --env mytd1 --env pd --env j /home/reevesp/bin/parallel --jobs 96 --env mytd1 --env pd --env j mytd1 >> "$pd"/"$j".distallel.txt;
    
    echo "$j: counting all alleles";
    >"$pd"/"$j".totallel.txt; #initialize output file with a header for count of all alleles within populations
    echo "$a" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env mytd2 --env pd --env k /home/reevesp/bin/parallel --jobs 96 --env mytd2 --env pd --env k mytd2 >> "$pd"/"$j".totallel.txt;
  done;


#count private alleles
grep ^'@' log.txt | tr '-' '_' | cut -d$'\t' -f1,3 | sort -t_ -k2,2n > freqs.txt;
grep ^'@' log.txt | tr '-' '_' | cut -d$'\t' -f1,2 | sort -t_ -k2,2n > counts.txt;

c=$(rev counts.txt | cut -d. -f2- | rev | uniq); #list of contig:site-ranges
>names.txt;
echo "$c" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env pd --env mypa /home/reevesp/bin/parallel --jobs 96 --env pd --env mypa mypa >> names.txt; #counting sub

#pad the file containing counts of private alleles with "?" for contig:site-range_readgroup combinations that weren't found
for bb in $b;
do >"$bb".privallel.txt; #create an output file for the readgroup
  dd="$bb"; export dd; #transfer iterator to variable that can be exported for gnu parallel
  echo "$a" | sed 's/#/@/g' | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env pd --env dd --env mypadpa \
                              /home/reevesp/bin/parallel --jobs 96 --env pd --env dd --env mypadpa mypadpa >> "$bb".privallel.txt;
done;

#sort and tab delimit output
for i in $b;
  do sort -t'_' -k1,1 "$i".distallel.txt | sort -t'_' -k2,2n | sed 's/^#//' | sed "1i position $i.distallel" | tr ' ' '\t' > "$i".2.distallel.txt;
    sort -t'_' -k1,1 "$i".totallel.txt | sort -t'_' -k2,2n | sed 's/^#//' | sed "1i position $i.totallel" | tr ' ' '\t'  > "$i".2.totallel.txt;
    sort -t'_' -k1,1 "$i".privallel.txt | sort -t'_' -k2,2n | sed 's/^@//' | sed "1i position $i.privallel" | tr ' ' '\t'  > "$i".2.privallel.txt;
  done;
#swap back to original filename
for i in $b; 
  do mv "$i".2.distallel.txt "$i".distallel.txt;
    mv "$i".2.totallel.txt "$i".totallel.txt;
    mv "$i".2.privallel.txt "$i".privallel.txt;
  done;
### END BASH ###


#consolidate files into 1, this part depends on number of read groups, is not universal
#you can tailor the columns cut depending on the number of read groups included
#for <= 6 read groups, it turns out the same
paste -d$'\t' global.distallel.txt RG_Z_*.distallel.txt global.totallel.txt RG_Z_*.totallel.txt global.privallel.txt RG_Z_*.privallel.txt \
  | cut -d$'\t' -f1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42 | sed 's/\.global//'> summary.txt;

#Pools differed in total number of reads. Scale proportionally by total read count relative to pool 51,
#which had the fewest
			L1
			pool:cols:scalar
			50:3,10,17:1.71
			51:4,11,18:1.00
			52:5,12,19:1.23
			53:6,13,20:1.83
			54:7,14,21:1.18
			55:8,15,22:1.20
			L1global
			2,9,16
			$3..$8=$2
			$10..$15=$9
			$17..$22=$16

cd /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1/Hs1pro1l1Finalhapx;

summ=$(tail -n +2 summary.txt | tr "\t" " "); #acquire summary.txt as a variable
h=$(head -1 summary.txt);
#adjust values for each read group using scaling values above
l1="50:3,10,17:1.71 51:4,11,18:1.00 52:5,12,19:1.23 53:6,13,20:1.83 54:7,14,21:1.18 55:8,15,22:1.20";
for i in $l1;
  do echo "$i";
    c=$(echo "$i" | cut -d: -f2); #columns to scale
    s=$(echo "$i" | cut -d: -f3); #scalar, divide by this
    #cycle through list of columns to scale
    for cc in $(echo "$c" | sed 's/,/ /g');
      do summ=$(echo "$summ" | awk -F$' ' -v cc=$cc -v s=$s '{$cc=$cc/s; print $0}');
      done;
  done;

#adjust the so-called 'global' value, which is just the sum of the read group values
l1global="2:3:8 9:10:15 16:17:22"; #columns to sum
for i in $l1global;
  do s=$(echo "$i" | cut -d: -f1); #column to put sum in
    st=$(echo "$i" | cut -d: -f2); #first column to sum
    se=$(echo "$i" | cut -d: -f3); #last column to sum
    summ=$(echo "$summ" | awk -F' ' -v s=$s -v st=$st -v se=$se '{for(j=st;j<=se;j++)x+=$j;$s=x; print $0; x=0}');
  done;

#write
echo "$h" > summaryscaled.txt;
echo "$summ" | awk -F' ' '{for(i=1; i<=NF; i++) {if($i=="0") $i="?"}; print $0}' | tr ' ' '\t' >> summaryscaled.txt;

### END OPTIONAL: COUNT TOTAL, DISTINCT, AND PRIVATE ALLELES ###



### IDENTIFY MICROHAPLOBLOCKS ACROSS HS1PRO1 ###

#Use hapxm.sh to find microhaploblocks in the bam files for each EL10 sco output by hapx -mb.
#By using -db and -va you can produce an
#optional file, mhendsvar.txt, that contains a variant optimized tiling path across the locus thru the
#microhaploblocks, which can be used directly as input into -u, to specify the universal microhaploblock set.
#hapxm option -va calculates a tiling path that collects the most variable (then, for ties, 
#most depth, then longest) microhaploblocks at each site.

#Hs1pro1l1, ~4 minutes on all nodes
cd /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1;
mkdir hapxmHs1Pro1L1;
cd hapxmHs1Pro1L1;
time hapxm.sh -b /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1/bam/Hs1pro1l1.finalaln.merge.bam \
            -ssh /home/reevesp/machines \
            -o hxm1 -db -va -s <(for i in $(seq 9889 1 15386); do echo 51jcf7180007742276:"$i"-"$i"; done;)

#use hapxm.sh to find microhaploblocks in the finalaln.bam.TMP file from each read group using the universal
#microhaplotypes discovered above. ~12 minutes on all nodes
time for i in $(seq 50 1 55);
  do echo $i;
    hapxm.sh -b /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1/bam/"$i"Hs1pro1l1.finalaln.bam.TMP \
            -ssh /home/reevesp/machines \
            -o hxm"$i"onhxm1 -u /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1/hapxmHs1Pro1L1/hxm1/mhendsvar.txt \
            -s <(for i in $(seq 9889 1 15386); do echo 51jcf7180007742276:"$i"-"$i"; done;)
  done;
### END IDENTIFY MICROHAPLOBLOCKS ACROSS HS1PRO1 ###




### FIND MICROHAPLOBLOCKS WHERE MAJOR ALLELE DIFFERS BETWEEN POOLS ###

cd /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1/hapxmHs1Pro1L1;
rhead=$(grep ^"#" hxm50onhxm1/hapxmlog.txt | cut -d' ' -f1-4);
#head=$(head -1 <<<"$rhead");
rhead=$(tail -n +2 <<<"$rhead");#remove col head from row heads

outmajall="";
outmajfreq="";
for i in $(seq 50 1 55);
  do echo $i;
    majall=$(grep ^"#" hxm"$i"onhxm1/hapxmlog.txt | tail -n +2 | cut -d' ' -f9 | cut -d: -f1); #extract the allele with highest frequency
    majfreq=$(grep ^"#" hxm"$i"onhxm1/hapxmlog.txt | tail -n +2 | cut -d' ' -f8 | cut -d: -f1); #extract the frequency of the major allele
    outmajall=$(paste -d: <(echo "$outmajall") <(echo "$majall"));
    outmajfreq=$(paste -d: <(echo "$outmajfreq") <(echo "$majfreq"));
  done;
outmajall=$(sed 's/^://' <<<"$outmajall"); #remove leading colon
outmajfreq=$(sed 's/^://' <<<"$outmajfreq"); #remove leading colon

#recode mhblocks from DNA sequence to integers
set -f; #you have to undo globbing before running this command since some values are solely asterisks
outrecode="";
while read l;
do m=$(echo "$l" | tr ':' '\n' | sed 's/^$/0/' | tr '\n' ':' | sed 's/:$//'); #substitute 0 when no allele was found
  a=$(echo "$m" | tr ':' '\n' | sort -u | grep -v 0); #list containing DNA sequences of alleles, excluding 0 (missing)
  j=1; #numeric allele name
  ll="$m";
  #if [[ "$a" == "" ]];
  #then outrecode+="$l"$'\n'; #deal with mhloci with no alleles by placing :::::
  #else
    for aa in $a; 
      do ll=${ll//"$aa"/"$j"}; #built-in bash replace obviates need to escape asterisks using sed
        j=$(( $j + 1 ));
      done;
     outrecode+="$ll"$'\n';
  #fi;
done <<<"$outmajall";
outrecode=$(sed '/^$/d' <<<"$outrecode"); #remove terminal blank line
set +f; #redo globbing

#assemble output
outp1=$(paste -d' ' <(echo "$rhead") <(echo "$outmajall") <(echo "$outrecode") <(echo "$outmajfreq")); #add row header to major allele list

#find all mh loci with different major alleles, label those rows with @
outp2="";
outp2=$(while read l;
do a=$(echo "$l" | cut -d' ' -f5 | tr ':' '\n'| sed '/^$/d' | sort | uniq | wc -l);
  if [[ "$a" > 1 ]];
  then echo "@$l";
  else echo "$l"
  fi;
done <<<"$outp1";
);

#write out the summary
echo "$outp2" > hxmsummary.txt;

### END FIND MICROHAPLOBLOCKS WHERE MAJOR ALLELE DIFFERS BETWEEN POOLS ###



### CREATE INPUT FILES FOR R ROUTINE MicrohaploblockHeatmap.r ###

cd /share/space/reevesp/patellifolia/AlleleMining/Hs1pro1/hapxmHs1Pro1L1;
c=$(cut -d' ' -f2-8 hxmsummary.txt); #get mhstart mhend mhlength seqs majalleles freqs for all mh majalleles
d=$(grep ^@ hxmsummary.txt | cut -d' ' -f2-8); #get mhstart mhend mhlength seqs majalleles freqs for mh majalleles that differ between pops

for k in "$c" "$d";
  do hdr="pop mhstart mhend length majallele freq seq";
    e="";
    while read l;
      do sqc=$(echo "$l" | cut -d' ' -f4); #allele sequence column
        mac=$(echo "$l" | cut -d' ' -f5); #major alleles column
        frc=$(echo "$l" | cut -d' ' -f6); #frequencies column
        j=1; #j indexes position in horizontal allele calls and freqs
        for i in $(seq 50 1 55);
          do rhead=$(echo -n "$i ";echo "$l" | cut -d' ' -f1-3); #row header labeled by population name
            sq=$(echo "$sqc" | cut -d: -f$j); #allele sequence for population $i
            ma=$(echo "$mac" | cut -d: -f$j); #major allele for population $i
            fr=$(echo "$frc" | cut -d: -f$j); #major allele frequency for population $i
            e+="$rhead $ma $fr $sq"$'\n';
    
            j=$(( $j + 1 ));
          done;
      done<<<"$k";
    e=$(sed '/^$/d' <<<"$e");
    e="$hdr"$'\n'"$e";
    
    
    if [[ "$k" == "$c" ]];
    then echo "$e" | tr ' ' '\t' > hxmsummaryRinAll.txt;
    elif [[ "$k" == "$d" ]]; 
    then echo "$e" | tr ' ' '\t' > hxmsummaryRinDiff.txt;
    fi;
  done;

### END CREATE INPUT FILES FOR R ROUTINE MicrohaploblockHeatmap.r ###



### PLOT MAJOR ALLELE MICROHAPLOBLOCK DIFFERENCES ###

#rsync everything back to local machine for work in R
cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/Hs1pro1;

### BEGIN R MicrohaploblockHeatmap.r ###
#install.packages("RColorBrewer")

options(error = recover)
rm(list=ls()) 
library("RColorBrewer")

# Functions #
	opacitybyallelefreq <- function(i,b,g)
	{
	  rr=rgb(b[1,i],b[2,i],b[3,i],alpha=g$freq[i]*255,max=255)
	  return(rr)
	}

	plotfeatures=function(f,fname,bbegin,eend)
	{
	  segments(f[1],eend+1,f[2],eend+1,col="black",lwd=20,lend=1)
	  segments(f[1],bbegin-1,f[1],eend+1,col="black",lwd=1,lend=1)
	  segments(f[2],bbegin-1,f[2],eend+1,col="black",lwd=1,lend=1)
	  text(f[1]+((f[2]-f[1])/2),eend+1,fname,col="white",cex=0.5,font=1)
	}
# End Functions #

# Main #
origpar=par() #gather initial settings
plottype=1 #1, fill empty space to next mh locus; 2, plot actual end points of mh

genelist=c("Hs1pro1L1")
for (gene in genelist)
{
	print(gene)
	
	if (gene=="Hs1pro1L1")
	{
		setwd("/Users/wichita/Desktop/telework/patellifolia/AlleleMining/Hs1pro1/hapxmHs1Pro1L1")
		pools=c("pat S", "pat T", "pat S", "pro H", "pro T", "web GC") #labels in order 50:55 (Hs1Pro1L1)
		labelpos=50:55
	} 	

	#import data table
	infilenames=c("hxmsummaryRinAll.txt","hxmsummaryRinDiff.txt")
	for (infilename in infilenames)
	{
		print(infilename)
		
		outfilename=paste(gene,gsub("hxmsummaryRin","",infilename),sep="")
		outfilename=gsub(".txt",".pdf",outfilename)

		g <- read.table(infilename, header=TRUE, sep="\t")

		#add a column that defines the start of the next microhaploblock allele for the current one
		nd=length(unique(g$pop)) #number of lines to delete from mhstart to shift it for mhnext
		shft=g$mhstart[(nd+1):length(g$mhstart)] #remove first nd elements from mhstart
		shft=c(shft,tail(g$mhend,nd))
		g$mhnext=shft #add column that will extend the stop point to the beginning of the next locus, for segment length

		#create a list of hex rgb colors with opacity proportional to allele frequency
		#rcbcolors=col2rgb(brewer.pal(n = 8, name = "Dark2")) #get hex codes for 8 RColorBrewer colors, convert to rgb
		#rcbcolors=cbind(col2rgb("black"), rcbcolors) #add black in position 1
		#b=rcbcolors[,g$majallele+1] #get the RColorBrewer Dark2 color for all major alleles, as defined by integer
		b=col2rgb(g$majallele+1) #color by major allele state, add 1 so that missing data (0) means black (r color value 1)
								 #this only works until there are more than 7 colors (r has only 8 numerically coded numbers)
		tpc=c() #initialize a vector to contain opacity adjusted color to plot for allele frequency
		for (i in 1:nrow(g))
		{
		  tpc=c(tpc,opacitybyallelefreq(i,b,g))
		}

		#add row to table containing opacity adjusted color
		g$freqcol=tpc

		#plot an empty graph and add appropriately colored and sized line segments
		par(mai=c(2.25,1.6,0.8,0.4)) #origpar$mai=c(1.0,1.6,0.8,0.4)

		pdf(outfilename, width=7, height=4.35)

		plot(g$mhstart,g$pop,type="n",yaxt="n",xlim=c(min(g$mhstart),max(g$mhend)),
			ylim=c(min(g$pop)-1,max(g$pop)+1), xlab="", ylab="") #set up plot with no points showing
		axis(2, at=labelpos, labels=pools, las=1) 

		if (plottype==1)
		{
			#display segments with augmented ends, to fill empty space to next locus
			segments(g$mhstart,g$pop,g$mhnext,g$pop,col=g$freqcol,lwd=30,lend=1) #add colored lines of appropriate length for locus
		} else
		if (plottype==2)
		{
			#display microhaploblock segments with length corresponding to actual end points
			segments(g$mhstart,g$pop,g$mhend,g$pop,col=g$freqcol,lwd=30,lend=1) #add colored lines of appropriate length for locus
		}

		#annotate graph with functional regions
		if (gene=="Hs1pro1L1")
		{
			fname="nematode feeding\nresponse" #name of feature
			f=c(12349,12951) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="general\nupregulation" #name of feature
			f=c(11392,11886) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="coding sequence" #name of feature
			f=c(12980,13828) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])
		}
		dev.off()
	} #infilenames
} #genelist
	
par(mai=origpar$mai)
# End Main #

### END R MicrohaploblockHeatmap.r ###

### END PLOT MAJOR ALLELE MICROHAPLOBLOCK DIFFERENCES ###
### END HS1PRO1 ALLELE MINING ###




### BVBTC1 ALLELE MINING ###
### DETERMINE CONTIGS CONTAINING BvBTC1 IN PATELLIFOLIA POOL ASSEMBLIES ###

module load blast+/2.9.0;

cd /home/pat.reeves/patellifolia/FlashedReadArchive;

#BvBTC1
for i in $(seq 50 1 55);
  do blastn -db /home/pat.reeves/patellifolia/FlashedReadArchive/"$i"fraFinal/"$i"Hs1pro1REVdref.fasta -query /home/pat.reeves/patellifolia/seq/BvBTC1mRNA.fa -out "$i"mapxBvBTC1.txt;
  done;

#consolidated output from blast
#all searches above return 1 hit in procumbens/webbiana, two hits in all three patellaris because it is tetraploid
			Query= HQ709091.1:2001-2131,2259-2297,2693-3192,4100-4261,4421-4554,
			4851-5012,10105-10285,10367-10802,10877-11699,11991-12338 Beta
			vulgaris subsp. vulgaris genotype KWS2320 bolting time control 1
			(BTC1) gene, complete cds

			Length=2916
																				  Score        E
			pool    Sequences producing significant alignments:                          (Bits)     Value
			50      jcf7180008080524                                              1070       0.0  
			50      jcf7180008061124                                              1048       0.0  
			51      jcf7180007719437                                              1070       0.0  
			51      jcf7180007731104                                              1048       0.0  
			52      jcf7180008029119                                              1070       0.0  
			52      jcf7180007952581                                              1048       0.0  
			53      jcf7180008837911                                              1046       0.0  
			54      jcf7180009178892                                              1046       0.0  
			55      jcf7180008533846                                              1046       0.0  


#Manual inspection of aligned sequences using sequencher showed shared indels can be used to
#distinguish orthologs.
#Within tetraploid patellaris one of the two sequences returned by blast is orthologous to
#the single sequence found in procumbens/webbiana

cd /home/pat.reeves/patellifolia/FlashedReadArchive;
sgrep '>jcf7180008080524_1_' 50fraFinal/50frasorted.fa | tr ' ' '\n' > 50jcf7180008080524.fa;
sgrep '>jcf7180008061124_1_' 50fraFinal/50frasorted.fa | tr ' ' '\n' > 50jcf7180008061124.fa;
sgrep '>jcf7180007719437_1_' 51fraFinal/51frasorted.fa | tr ' ' '\n' > 51jcf7180007719437.fa;
sgrep '>jcf7180007731104_1_' 51fraFinal/51frasorted.fa | tr ' ' '\n' > 51jcf7180007731104.fa;
sgrep '>jcf7180008029119_1_' 52fraFinal/52frasorted.fa | tr ' ' '\n' > 52jcf7180008029119.fa;
sgrep '>jcf7180007952581_1_' 52fraFinal/52frasorted.fa | tr ' ' '\n' > 52jcf7180007952581.fa;
sgrep '>jcf7180008837911_1_' 53fraFinal/53frasorted.fa | tr ' ' '\n' > 53jcf7180008837911.fa;
sgrep '>jcf7180009178892_1_' 54fraFinal/54frasorted.fa | tr ' ' '\n' > 54jcf7180009178892.fa;
sgrep '>jcf7180008533846_1_' 55fraFinal/55frasorted.fa | tr ' ' '\n' > 55jcf7180008533846.fa;
mv 5[0-5]jcf*.fa /home/pat.reeves/patellifolia/AlleleMining/BTC1;




			BTC1 locus 1 orthologs, and spans containing the BTC1 gene + 2kb upstream are:
			50jcf7180008061124_5108-17079
			51jcf7180007731104_4958-16929
			*52jcf7180007952581_15244-27196
			53jcf7180008837911_32968-44910
			54jcf7180009178892_5755-17712
			55jcf7180008533846_371-12324
			
			* = proper sense orientation, use this one as single reference gene for mapping against

			BTC1 locus 2 orthologs, and spans containing BTC1L2 gene (not used here) are:
			50jcf7180008080524
			51jcf7180007719437
			52jcf7180008029119
						
#After identifying which contigs contain BTC1 locus 1 in each pool assembly, and the
#coordinates within those contigs that contain the BTC1 gene, extract the phased reads
#that map to that region from the FRA map. Just use samtools for this since the reads are
#already phased and quality filtered by virtue of being in the map
			
cd /home/pat.reeves/patellifolia/AlleleMining/BTC1;
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/50fraFinal/50fra.bam "jcf7180008061124:5108-17079" > 50BTC1L1.bam; #11972bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/51fraFinal/51fra.bam "jcf7180007731104:4958-16929" > 51BTC1L1.bam; #11972bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/52fraFinal/52fra.bam "jcf7180007952581:15244-27196" > 52BTC1L1.bam; #11953bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/53fraFinal/53fra.bam "jcf7180008837911:32968-44910" > 53BTC1L1.bam; #11943bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/54fraFinal/54fra.bam "jcf7180009178892:5755-17712" > 54BTC1L1.bam; #11958bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/55fraFinal/55fra.bam "jcf7180008533846:371-12324" > 55BTC1L1.bam; #11954bp

#extract phased reads from bam files into a single fasta file
seq 50 1 55 | parallel --keep-order 'echo {}; samtools fasta {}BTC1L1.bam > {}BTC1L1.fa';

			50
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 2535 reads
			51
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 1494 reads
			52
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 2072 reads
			53
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 5385 reads
			54
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 3443 reads
			55
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 3609 reads
			

#manual depth calculation from samtools stats 5[0-5]BTC1l1.finalaln.bam.TMP | grep 'sequences\|length' | grep -v ^#;
			pool	sequences	longest frag	total length	phased read count	phased read total length	cistron length	avg phased read depth
			50		2535		35291			737948			2534				702657						11972			58.69
			51		1494		104155			515814			1493				411659						11972			34.39
			52		2072		32297			606242			2071				573945						11953			48.02
			53		5385		58036			1520104			5384				1462068						11943			122.42
			54		3443		25987			976417			3442				950430						11958			79.48
			55		3609		25414			942496			3608				917082						11954			76.72


#map these phased reads to a single reference gene sequence (choose 52jcf7180007952581,
#the only one in sense orientation)
#for BTC1L1 make 52jcf7180007952581_ref.txt from 52jcf7180007952581.fa, indexed with bwa index.
#you must rename the single reference sequence from "jcf7180007952581_1_32297" because a
#sequence by that name is included in the fasta files from above. It will cause hapx to
#fail in the calculation of $tmpf
cd /home/pat.reeves/patellifolia/AlleleMining/BTC1;
mkdir bam ref;
sed 's/jcf7180007952581_1_32297/52jcf7180007952581_ref/' 52jcf7180007952581.fa > ref/52jcf7180007952581_ref.txt;
bwa index ref/52jcf7180007952581_ref.txt;
pd=$(pwd);
t=40; #threads
for j in 50BTC1L1.fa 51BTC1L1.fa 52BTC1L1.fa 53BTC1L1.fa 54BTC1L1.fa 55BTC1L1.fa;
    do rg=$(cut -c1-2 <<<"$j"); #get the read group name 5[0-5]
      echo "$rg";
      readgroup="@RG\tID:$rg\tSM:$rg\tPL:illumina\tLB:na\tPU:na";
      bwa mem -R "$readgroup" -t "$t" "$pd"/ref/52jcf7180007952581_ref.txt "$pd"/"$j" | \
                  samtools sort -O BAM --threads "$t" | \
                  samtools view -F 2048 -O BAM > "$pd"/bam/"$rg"BTC1l1.finalaln.bam.TMP; #-F 2048 excludes supplementary alignments
                  samtools index "$pd"/bam/"$rg"BTC1l1.finalaln.bam.TMP;
    done;
sambamba merge -t40 -p "$pd"/bam/BTC1l1.finalaln.merge.bam "$pd"/bam/5[0-5]BTC1l1.finalaln.bam.TMP;


### END DETERMINE CONTIGS CONTAINING BvBTC1 IN PATELLIFOLIA POOL ASSEMBLIES ###




### OPTIONAL: COUNT TOTAL, DISTINCT, AND PRIVATE ALLELES ###
#rsync results from HPC ceres to HPC blip so you can get access to multi-node gnu parallel

#run a final hapx to count microhaplotypes at each position (takes ~4.5 hrs)
cd /share/space/reevesp/patellifolia/AlleleMining/BTC1; #go here just to get the working directory specified
pd=$(pwd);
time hapx.sh -r "$pd"/ref/52jcf7180007952581_ref.txt \
        -b "$pd"/bam/BTC1l1.finalaln.merge.bam \
        -o BTC1l1Finalhapx \
        -f 0 -q 1 -x -ssh /home/reevesp/machines \
        -s <(for i in $(seq 14000 1 28000); do echo 52jcf7180007952581_ref:"$i"-"$i"; done;); #set window to +- ~1000bp around locus in ref

#postprocess haploblock counts in log.txt into something plottable
cd /share/space/reevesp/patellifolia/AlleleMining/BTC1/BTC1l1Finalhapx;

### BEGIN BASH ###
mytd1() {
        i=$1; #a position on the reference is incoming
        a=$(grep "$i"\\."$k" "$pd"/numblocks.txt | cut -d: -f4); #number of unique haploblock read pairs mapped to position i in the readgroup
        if [[ "$a" == "" ]]; then a="?"; fi; #if no data at position i set number of haploblock read pairs to ?
        echo "$i"."$k $a"; #report result to parallel statement
}
export -f mytd1;

mytd2() {
        i=$1; #a position on the reference is incoming
        b=$(grep "$i"\\."$k" "$pd"/numblocks.txt | cut -d$'\t' -f2 | cut -d: -f1); #total number haploblock read pairs mapped to position i in the readgroup
        if [[ "$b" == "" ]]; then b="?"; fi; #if no data at position i set to ?
        echo "$i"."$k $b"; #report result to parallel statement
}
export -f mytd2;

mypa() {
       i=$1; #contig:site-range
       j=$(grep "$i" "$pd"/counts.txt | grep -v global);
       names=$(cut -d$'\t' -f1 <<<"$j" | sed 's/$/ 0/');
       counts=$(cut -d$'\t' -f2 <<<"$j");
       nc=$(head -1 <<<"$counts" | awk -F: '{print NF}'); #number of alleles
       nr=$(wc -l <<<"$counts"); #number of read groups
       nz=$(( $nr - 1 )); #number of zeroes needed in a column of allele counts for a private allele to exist
       for m in $(seq 1 1 $nc);
         do col=$(cut -d: -f$m <<<"$counts"); #extract the column of allele counts
           fl=$(grep -n -v 0 <<<"$col"); #show line numbers where alleles are present
           
           #if only one line of the column of data has alleles, it is a private allele
           if [[ $(echo "$fl" | wc -l) == 1 ]];
           then r=$(cut -d: -f1 <<<"$fl"); #row number of private allele
             names=$(awk -F' ' -v r=$r '{if (NR==r) $2++; print}' <<<"$names"); #index up by one the row with the private allele
           fi;
         done;
         
       #report result to parallel statement
       echo "$names";
}
export -f mypa;

mypadpa() {
          aa=$1;
          cc=$(grep "$aa.$dd" "$pd"/names.txt);
          if [[ "$cc" == "" ]];
          then echo "$aa.$dd ?"; #if readgroup not found print a ?
          else echo "$cc";
          fi;
}
export -f mypadpa;


pd=$(pwd); export pd;

#count total and distinct haploblocks
grep ^# log.txt | tr '-' '_' > numblocks.txt;
a=$(rev numblocks.txt | cut -d. -f2- | rev | uniq); #list of contig:site-ranges
b=$(rev numblocks.txt | cut -d. -f1 | rev | cut -d$'\t' -f1 | sort -u); export b; #determine set of possible read groups

#iterate over read groups to count haploblocks
for j in $b;
  do k="$j"; export k; #do this so iterator can be sent to nodes using --sshloginfile
    echo "$j: counting distinct alleles";
    >"$pd"/"$j".distallel.txt; #initialize output file with a header for count of distinct alleles within populations
    echo "$a" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env mytd1 --env pd --env k /home/reevesp/bin/parallel --jobs 96 --env mytd1 --env pd --env k mytd1 >> "$pd"/"$j".distallel.txt;
    #echo "$a" | parallel --jobs 1 --pipe -N960 --env mytd1 --env pd --env j /home/reevesp/bin/parallel --jobs 96 --env mytd1 --env pd --env j mytd1 >> "$pd"/"$j".distallel.txt;
    
    echo "$j: counting all alleles";
    >"$pd"/"$j".totallel.txt; #initialize output file with a header for count of all alleles within populations
    echo "$a" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env mytd2 --env pd --env k /home/reevesp/bin/parallel --jobs 96 --env mytd2 --env pd --env k mytd2 >> "$pd"/"$j".totallel.txt;
  done;


#count private alleles
grep ^'@' log.txt | tr '-' '_' | cut -d$'\t' -f1,3 | sort -t_ -k2,2n > freqs.txt;
grep ^'@' log.txt | tr '-' '_' | cut -d$'\t' -f1,2 | sort -t_ -k2,2n > counts.txt;

c=$(rev counts.txt | cut -d. -f2- | rev | uniq); #list of contig:site-ranges
>names.txt;
echo "$c" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env pd --env mypa /home/reevesp/bin/parallel --jobs 96 --env pd --env mypa mypa >> names.txt; #counting sub

#pad the file containing counts of private alleles with "?" for contig:site-range_readgroup combinations that weren't found
for bb in $b;
do >"$bb".privallel.txt; #create an output file for the readgroup
  dd="$bb"; export dd; #transfer iterator to variable that can be exported for gnu parallel
  echo "$a" | sed 's/#/@/g' | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env pd --env dd --env mypadpa \
                              /home/reevesp/bin/parallel --jobs 96 --env pd --env dd --env mypadpa mypadpa >> "$bb".privallel.txt;
done;

#sort and tab delimit output
for i in $b;
  do sort -t'_' -k1,1 "$i".distallel.txt | sort -t'_' -k2,2n | sed 's/^#//' | sed "1i position $i.distallel" | tr ' ' '\t' > "$i".2.distallel.txt;
    sort -t'_' -k1,1 "$i".totallel.txt | sort -t'_' -k2,2n | sed 's/^#//' | sed "1i position $i.totallel" | tr ' ' '\t'  > "$i".2.totallel.txt;
    sort -t'_' -k1,1 "$i".privallel.txt | sort -t'_' -k2,2n | sed 's/^@//' | sed "1i position $i.privallel" | tr ' ' '\t'  > "$i".2.privallel.txt;
  done;
#swap back to original filename
for i in $b; 
  do mv "$i".2.distallel.txt "$i".distallel.txt;
    mv "$i".2.totallel.txt "$i".totallel.txt;
    mv "$i".2.privallel.txt "$i".privallel.txt;
  done;
### END BASH ###


#consolidate files into 1, this part depends on number of read groups, is not universal
#you can tailor the columns cut depending on the number of read groups included
#for <= 6 read groups, it turns out the same
paste -d$'\t' global.distallel.txt RG_Z_*.distallel.txt global.totallel.txt RG_Z_*.totallel.txt global.privallel.txt RG_Z_*.privallel.txt \
  | cut -d$'\t' -f1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42 | sed 's/\.global//'> summary.txt;

#Pools differed in total number of reads. Scale proportionally by total read count relative to pool 51,
#which had the fewest
			L1
			pool:cols:scalar
			50:3,10,17:1.71
			51:4,11,18:1.00
			52:5,12,19:1.23
			53:6,13,20:1.83
			54:7,14,21:1.18
			55:8,15,22:1.20
			L1global
			2,9,16
			$3..$8=$2
			$10..$15=$9
			$17..$22=$16

cd /share/space/reevesp/patellifolia/AlleleMining/BTC1/BTC1l1Finalhapx;

summ=$(tail -n +2 summary.txt | tr "\t" " "); #acquire summary.txt as a variable
h=$(head -1 summary.txt);
#adjust values for each read group using scaling values above
l1="50:3,10,17:1.71 51:4,11,18:1.00 52:5,12,19:1.23 53:6,13,20:1.83 54:7,14,21:1.18 55:8,15,22:1.20";
for i in $l1;
  do echo "$i";
    c=$(echo "$i" | cut -d: -f2); #columns to scale
    s=$(echo "$i" | cut -d: -f3); #scalar, divide by this
    #cycle through list of columns to scale
    for cc in $(echo "$c" | sed 's/,/ /g');
      do summ=$(echo "$summ" | awk -F$' ' -v cc=$cc -v s=$s '{$cc=$cc/s; print $0}');
      done;
  done;

#adjust the so-called 'global' value, which is just the sum of the read group values
l1global="2:3:8 9:10:15 16:17:22"; #columns to sum
for i in $l1global;
  do s=$(echo "$i" | cut -d: -f1); #column to put sum in
    st=$(echo "$i" | cut -d: -f2); #first column to sum
    se=$(echo "$i" | cut -d: -f3); #last column to sum
    summ=$(echo "$summ" | awk -F' ' -v s=$s -v st=$st -v se=$se '{for(j=st;j<=se;j++)x+=$j;$s=x; print $0; x=0}');
  done;

#write
echo "$h" > summaryscaled.txt;
echo "$summ" | awk -F' ' '{for(i=1; i<=NF; i++) {if($i=="0") $i="?"}; print $0}' | tr ' ' '\t' >> summaryscaled.txt;

### END OPTIONAL: COUNT TOTAL, DISTINCT, AND PRIVATE ALLELES ###



### IDENTIFY MICROHAPLOBLOCKS ACROSS BvBTC1 ###

#Use hapxm.sh to find microhaploblocks in the bam files for each EL10 sco output by hapx -mb.
#By using -db and -va you can produce an
#optional file, mhendsvar.txt, that contains a variant optimized tiling path across the locus thru the
#microhaploblocks, which can be used directly as input into -u, to specify the universal microhaploblock set.
#hapxm option -va calculates a tiling path that collects the most variable (then, for ties, 
#most depth, then longest) microhaploblocks at each site.

#BTC1l1, ~11 minutes on all nodes
#determine range to consider
cd /share/space/reevesp/patellifolia/AlleleMining/BTC1/bam;
samtools view BTC1l1.finalaln.merge.bam | grep -v ^'jcf' | awk -F$'\t' '$4!=0{print $0}' | cut -d$'\t' -f4 | sort -n | head; #leftmost mapping position = 14761
samtools view BTC1l1.finalaln.merge.bam | grep -v ^'jcf' | awk -F$'\t' '$4!=0{print $0}' | cut -d$'\t' -f4 | sort -n | tail; #rightmost mapping position = 27253 + 150 = 27403

cd /share/space/reevesp/patellifolia/AlleleMining/BTC1;
mkdir hapxmBTC1L1;
cd hapxmBTC1L1;
time hapxm.sh -b /share/space/reevesp/patellifolia/AlleleMining/BTC1/bam/BTC1l1.finalaln.merge.bam \
            -ssh /home/reevesp/machines \
            -o hxm1 -db -va -s <(for i in $(seq 14761 1 27403); do echo 52jcf7180007952581_ref:"$i"-"$i"; done;)

#use hapxm.sh to find microhaploblocks in the finalaln.bam.TMP file from each read group using the universal
#microhaplotypes discovered above. ~23 minutes on all nodes
time for i in $(seq 50 1 55);
  do echo $i;
    hapxm.sh -b /share/space/reevesp/patellifolia/AlleleMining/BTC1/bam/"$i"BTC1l1.finalaln.bam.TMP \
            -ssh /home/reevesp/machines \
            -o hxm"$i"onhxm1 -u /share/space/reevesp/patellifolia/AlleleMining/BTC1/hapxmBTC1L1/hxm1/mhendsvar.txt \
            -s <(for i in $(seq 14761 1 27403); do echo 52jcf7180007952581_ref:"$i"-"$i"; done;)
  done;
### END IDENTIFY MICROHAPLOBLOCKS ACROSS BvBTC1 ###



### FIND MICROHAPLOBLOCKS WHERE MAJOR ALLELE DIFFERS BETWEEN POOLS ###

cd /share/space/reevesp/patellifolia/AlleleMining/BTC1/hapxmBTC1L1;
rhead=$(grep ^"#" hxm50onhxm1/hapxmlog.txt | cut -d' ' -f1-4);
rhead=$(tail -n +2 <<<"$rhead");#remove col head from row heads

outmajall="";
outmajfreq="";
for i in $(seq 50 1 55);
  do echo $i;
    majall=$(grep ^"#" hxm"$i"onhxm1/hapxmlog.txt | tail -n +2 | cut -d' ' -f9 | cut -d: -f1); #extract the allele with highest frequency
    majfreq=$(grep ^"#" hxm"$i"onhxm1/hapxmlog.txt | tail -n +2 | cut -d' ' -f8 | cut -d: -f1); #extract the frequency of the major allele
    outmajall=$(paste -d: <(echo "$outmajall") <(echo "$majall"));
    outmajfreq=$(paste -d: <(echo "$outmajfreq") <(echo "$majfreq"));
  done;
outmajall=$(sed 's/^://' <<<"$outmajall"); #remove leading colon
outmajfreq=$(sed 's/^://' <<<"$outmajfreq"); #remove leading colon

#recode mhblocks from DNA sequence to integers
set -f; #you have to undo globbing before running this command since some values are solely asterisks
outrecode="";
while read l;
do m=$(echo "$l" | tr ':' '\n' | sed 's/^$/0/' | tr '\n' ':' | sed 's/:$//'); #substitute 0 when no allele was found
  a=$(echo "$m" | tr ':' '\n' | sort -u | grep -v 0); #list containing DNA sequences of alleles, excluding 0 (missing)
  j=1; #numeric allele name
  ll="$m";
  #if [[ "$a" == "" ]];
  #then outrecode+="$l"$'\n'; #deal with mhloci with no alleles by placing :::::
  #else
    for aa in $a; 
      do ll=${ll//"$aa"/"$j"}; #built-in bash replace obviates need to escape asterisks using sed
        j=$(( $j + 1 ));
      done;
     outrecode+="$ll"$'\n';
  #fi;
done <<<"$outmajall";
outrecode=$(sed '/^$/d' <<<"$outrecode"); #remove terminal blank line
set +f; #redo globbing

#assemble output
outp1=$(paste -d' ' <(echo "$rhead") <(echo "$outmajall") <(echo "$outrecode") <(echo "$outmajfreq")); #add row header to major allele list

#find all mh loci with different major alleles, label those rows with @
outp2="";
outp2=$(while read l;
do a=$(echo "$l" | cut -d' ' -f5 | tr ':' '\n'| sed '/^$/d' | sort | uniq | wc -l);
  if [[ "$a" > 1 ]];
  then echo "@$l";
  else echo "$l"
  fi;
done <<<"$outp1";
);

#write out the summary
echo "$outp2" > hxmsummary.txt;

### END FIND MICROHAPLOBLOCKS WHERE MAJOR ALLELE DIFFERS BETWEEN POOLS ###



### CREATE INPUT FILES FOR R ROUTINE MicrohaploblockHeatmap.r ###

cd /share/space/reevesp/patellifolia/AlleleMining/BTC1/hapxmBTC1L1;
c=$(cut -d' ' -f2-8 hxmsummary.txt); #get mhstart mhend mhlength seqs majalleles freqs for all mh majalleles
d=$(grep ^@ hxmsummary.txt | cut -d' ' -f2-8); #get mhstart mhend mhlength seqs majalleles freqs for mh majalleles that differ between pops

for k in "$c" "$d";
  do hdr="pop mhstart mhend length majallele freq seq";
    e="";
    while read l;
      do sqc=$(echo "$l" | cut -d' ' -f4); #allele sequence column
        mac=$(echo "$l" | cut -d' ' -f5); #major alleles column
        frc=$(echo "$l" | cut -d' ' -f6); #frequencies column
        j=1; #j indexes position in horizontal allele calls and freqs
        for i in $(seq 50 1 55);
          do rhead=$(echo -n "$i ";echo "$l" | cut -d' ' -f1-3); #row header labeled by population name
            sq=$(echo "$sqc" | cut -d: -f$j); #allele sequence for population $i
            ma=$(echo "$mac" | cut -d: -f$j); #major allele for population $i
            fr=$(echo "$frc" | cut -d: -f$j); #major allele frequency for population $i
            e+="$rhead $ma $fr $sq"$'\n';
    
            j=$(( $j + 1 ));
          done;
      done<<<"$k";
    e=$(sed '/^$/d' <<<"$e");
    e="$hdr"$'\n'"$e";
    
    
    if [[ "$k" == "$c" ]];
    then echo "$e" | tr ' ' '\t' > hxmsummaryRinAll.txt;
    elif [[ "$k" == "$d" ]]; 
    then echo "$e" | tr ' ' '\t' > hxmsummaryRinDiff.txt;
    fi;
  done;

### END CREATE INPUT FILES FOR R ROUTINE MicrohaploblockHeatmap.r ###



### CALCULATE DESCRIPTIVE STATISTICS FOR MICROHAPLOBLOCKS IN VARIANT-RICH TILING ARRAY ###

#myst() calculates basic stats
myst() {
  awk '{
    sum = 0;
    M = 0;
    S = 0;
    for (k=1; k <= NF; k++) {
      sum += $k;
      x = $k;
      oldM = M;
      M = M + ((x - M)/k);
      S = S + (x - M)*(x - oldM);
    }
    var = S/(NF - 1);
    print "n=" NF " mean=" sum/(NF) " var=" var " sd=" sqrt(var);
  }' <<< $1
}
export myst;

cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/BTC1/hapxmBTC1L1;

#Choose boundaries by visual inspection of hxmsummaryRinAll.txt and hxmsummaryRinDiff.txt in order to
#trim off long mhs and poor alignment regions from ends, which are artifacts. Choice is subjective.
for i in {50..55};
  do a=$(sed -n '/^#52jcf7180007952581_ref 15249 15250/,$p' hxm"$i"onhxm1/hapxmlog.txt | \
         sed -e '/^#52jcf7180007952581_ref 27163 27163/,$d'); #sed -n $p prints lines after including match, -e $d deletes lines after including match
     b=$(echo "$a" | awk -F' ' '$4!=1{print $0}'); #exclude loci with mhlength=1
     #b=$(echo "$b" | awk -F' ' '$6!=1{print $0}'); #exclude loci with only 1 variant (don't use this)
     nmh=$(echo "$b" | wc -l); #number of mhs across minimum tiling path
     
     #calc stats
     mhl=$(echo "$b" | cut -d' ' -f4 | tr "\n" " " | sed 's/ $//'); #stats on mhlength
     mhlengthst=$(myst "$mhl");
     ns=$(echo "$b" | cut -d' ' -f5 | tr "\n" " " | sed 's/ $//'); #stats on numseq
     numseqst=$(myst "$ns");
     na=$(echo "$b" | cut -d' ' -f6 | tr "\n" " " | sed 's/ $//'); #stats on numalleles
     numallelesst=$(myst "$na");
     
     #calc mean length of mhs
     #c=$(echo "$b" | awk -F' ' '{s+=$4} END {print s}'); #sum of mhlength
     #ml=$(echo "scale=4; $c / $nmh" | bc);
     
     echo "$i $mhlengthst | $numseqst | $numallelesst"
  done;

#result, basic stats, mh with 1 variant included:
			50 n=3220 mean=2.99379 var=1.92416 sd=1.38714 | n=3220 mean=35.3752 var=43.1354 sd=6.56775 | n=3220 mean=2.0736 var=0.282559 sd=0.531562
			51 n=3220 mean=2.99379 var=1.92416 sd=1.38714 | n=3220 mean=19.9317 var=24.6682 sd=4.96671 | n=3220 mean=1.64938 var=0.381842 sd=0.617933
			52 n=3220 mean=2.99379 var=1.92416 sd=1.38714 | n=3220 mean=30.1696 var=41.7066 sd=6.45806 | n=3220 mean=1.99689 var=0.389552 sd=0.624141
			53 n=3220 mean=2.99379 var=1.92416 sd=1.38714 | n=3220 mean=70.923 var=121.674 sd=11.0306 | n=3220 mean=2.92019 var=0.883656 sd=0.94003
			54 n=3220 mean=2.99379 var=1.92416 sd=1.38714 | n=3220 mean=47.5898 var=71.0398 sd=8.42851 | n=3220 mean=2.37298 var=0.477493 sd=0.691009
			55 n=3220 mean=2.99379 var=1.92416 sd=1.38714 | n=3220 mean=49.8099 var=82.2136 sd=9.06717 | n=3220 mean=2.23292 var=0.517338 sd=0.719262


#calculate the proportion of loci containing only indel variants
for i in {50..55};
  do 
     a=$(sed -n '/^#52jcf7180007952581_ref 15249 15250/,$p' hxm"$i"onhxm1/hapxmlog.txt | \
         sed -e '/^#52jcf7180007952581_ref 27163 27163/,$d'); #sed -n $p prints lines after including match, -e $d deletes lines after including match
     b=$(echo "$a" | awk -F' ' '$4!=1{print $0}'); #exclude lines with mhlength=1
     c=$(echo "$b" | awk -F' ' '$6!=1{print $0}' | cut -d' ' -f9); #collect variants, excluding loci with only 1 variant
     #c=$(echo "$b" | cut -d' ' -f9); #collect variants, including loci with only 1 variant
     
     #count loci that only contain indel variants
     indelloc=0; #number of loci that contain only indel type variants
     snploc=0
     while read d;
       do e=$(echo "$d" | tr ':' '\n' | awk NF=NF FS= | rs -c' ' -C' ' -T | sed 's/ //g' \
              | sed 's/\*//g' | grep -v ^"A\+"$ | grep -v ^"C\+"$ | grep -v ^"G\+"$ | grep -v ^"T\+"$ \
              | wc -l); #split to rows, add space btw each character, transpose, remove space delimiter
                        #delete asterisks, remove lines that are all repeated ACG or T
                        #count number of lines remaining (if >0 the locus contains at least 1 SNP variant
                        #if =0 the locus contains only indel variants)
          if (( $e == 0 )); then (( indelloc++ )); fi;
          if (( $e > 0 )); then (( snploc++ )); fi;
       done <<<"$c";
       
       totloc=$(( indelloc + snploc ));
       indelpro=$(echo "scale=4; $indelloc / $totloc" | bc);
       snppro=$(echo "scale=4; $snploc / $totloc" | bc);
     
     echo "$i $indelloc $indelpro $snploc $snppro $totloc";
  done;

#results, single variant loci excluded:
			50 2885 .9948 15 .0051 2900
			51 1837 .9919 15 .0080 1852
			52 2606 .9878 32 .0121 2638
			53 2977 .9376 198 .0623 3175
			54 2972 .9665 103 .0334 3075
			55 2764 .9674 93 .0325 2857

#repeat calculating the proportion of loci containing only indel variants for the combined dataset
cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/BTC1/hapxmBTC1L1;
for i in hxm1;
  do 
     a=$(sed -n '/^#52jcf7180007952581_ref 15249 15250/,$p' "$i"/hapxmlogvar.txt | \
         sed -e '/^#52jcf7180007952581_ref 27163 27163/,$d'); #sed -n $p prints lines after including match, -e $d deletes lines after including match
     b=$(echo "$a" | awk -F' ' '$4!=1{print $0}'); #exclude lines with mhlength=1
     c=$(echo "$b" | awk -F' ' '$6!=1{print $0}' | cut -d' ' -f9); #collect variants, excluding loci with only 1 variant
     #c=$(echo "$b" | cut -d' ' -f9); #collect variants, including loci with only 1 variant
     
     #count loci that only contain indel variants
     indelloc=0; #number of loci that contain only indel type variants
     snploc=0
     while read d;
       do e=$(echo "$d" | tr ':' '\n' | awk NF=NF FS= | rs -c' ' -C' ' -T | sed 's/ //g' \
              | sed 's/\*//g' | grep -v ^"A\+"$ | grep -v ^"C\+"$ | grep -v ^"G\+"$ | grep -v ^"T\+"$ \
              | wc -l); #split to rows, add space btw each character, transpose, remove space delimiter
                        #delete asterisks, remove lines that are all repeated ACG or T
                        #count number of lines remaining (if >0 the locus contains at least 1 SNP variant
                        #if =0 the locus contains only indel variants)
          if (( $e == 0 )); then (( indelloc++ )); fi;
          if (( $e > 0 )); then (( snploc++ )); fi;
       done <<<"$c";
       
       totloc=$(( indelloc + snploc ));
       indelpro=$(echo "scale=4; $indelloc / $totloc" | bc);
       snppro=$(echo "scale=4; $snploc / $totloc" | bc);
     
     echo "$i $indelloc $indelpro $snploc $snppro $totloc";
  done;


#results, single variant loci excluded:
			hxm1 2651 .8232 569 .1767 3220

#calculate major variant frequency
cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/BTC1/hapxmBTC1L1;
a=$(sed -n '/^#52jcf7180007952581_ref 15249 15250/,$p' hxmsummary.txt | \
    sed -e '/^#52jcf7180007952581_ref 27163 27163/,$d'); #sed -n $p prints lines after including match, -e $d deletes lines after including match
b=$(echo "$a" | awk -F' ' '$4!=1{print $0}'); #exclude lines with mhlength=1
c=$(echo "$b" | cut -d' ' -f7); #collect major variant frequencies, including loci with only 1 variant

#calc stats
for i in {50..55};
  do j=$(( $i - 49 )); #calculate a column index for each pool
    vf=$(echo "$c" | cut -d: -f"$j" | grep -v ^0$ | tr "\n" " " | sed 's/ $//'); #get variant freqs, remove zeroes (no variant called)
    varfreqs=$(myst "$vf");
    echo "$i $varfreqs";
  done;
			50 n=3220 mean=0.890574 var=0.00393507 sd=0.0627301
			51 n=3220 mean=0.909875 var=0.00904554 sd=0.095108
			52 n=3220 mean=0.894469 var=0.00524432 sd=0.0724177
			53 n=3220 mean=0.866678 var=0.00564243 sd=0.0751161
			54 n=3220 mean=0.875224 var=0.00526525 sd=0.0725621
			55 n=3214 mean=0.913399 var=0.00333679 sd=0.057765

### END CALCULATE DESCRIPTIVE STATISTICS FOR MICROHAPLOBLOCKS IN VARIANT-RICH TILING ARRAY ###



### PLOT MAJOR ALLELE MICROHAPLOBLOCK DIFFERENCES ###

#rsync everything back to local machine for work in R
cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/BTC1;

### BEGIN R MicrohaploblockHeatmap.r ###
#install.packages("RColorBrewer")

options(error = recover)
rm(list=ls()) 
library("RColorBrewer")

# Functions #
	opacitybyallelefreq <- function(i,b,g)
	{
	  rr=rgb(b[1,i],b[2,i],b[3,i],alpha=g$freq[i]*255,max=255)
	  return(rr)
	}

	plotfeatures=function(f,fname,bbegin,eend)
	{
	  segments(f[1],eend+1,f[2],eend+1,col="black",lwd=20,lend=1)
	  segments(f[1],bbegin-1,f[1],eend+1,col="black",lwd=1,lend=1)
	  segments(f[2],bbegin-1,f[2],eend+1,col="black",lwd=1,lend=1)
	  text(f[1]+((f[2]-f[1])/2),eend+1,fname,col="white",cex=0.5,font=1)
	}
# End Functions #

# Main #
origpar=par() #gather initial settings
plottype=1 #1, fill empty space to next mh locus; 2, plot actual end points of mh

genelist=c("BTC1L1")
for (gene in genelist)
{
	print(gene)
	
	if (gene=="BTC1L1")
	{
		setwd("/Users/wichita/Desktop/telework/patellifolia/AlleleMining/BTC1/hapxmBTC1L1")
		pools=c("pat S", "pat T", "pat S", "pro H", "pro T", "web GC") #labels in order 50:55 (BTC1L1)
		labelpos=50:55
	} 	

	#import data table
	infilenames=c("hxmsummaryRinAll.txt","hxmsummaryRinDiff.txt")
	for (infilename in infilenames)
	{
		print(infilename)
		
		outfilename=paste(gene,gsub("hxmsummaryRin","",infilename),sep="")
		outfilename=gsub(".txt",".pdf",outfilename)

		g <- read.table(infilename, header=TRUE, sep="\t")

		#add a column that defines the start of the next microhaploblock allele for the current one
		nd=length(unique(g$pop)) #number of lines to delete from mhstart to shift it for mhnext
		shft=g$mhstart[(nd+1):length(g$mhstart)] #remove first nd elements from mhstart
		shft=c(shft,tail(g$mhend,nd))
		g$mhnext=shft #add column that will extend the stop point to the beginning of the next locus, for segment length

		#create a list of hex rgb colors with opacity proportional to allele frequency
		#rcbcolors=col2rgb(brewer.pal(n = 8, name = "Dark2")) #get hex codes for 8 RColorBrewer colors, convert to rgb
		#rcbcolors=cbind(col2rgb("black"), rcbcolors) #add black in position 1
		#b=rcbcolors[,g$majallele+1] #get the RColorBrewer Dark2 color for all major alleles, as defined by integer
		b=col2rgb(g$majallele+1) #color by major allele state, add 1 so that missing data (0) means black (r color value 1)
								 #this only works until there are more than 7 colors (r has only 8 numerically coded numbers)
		tpc=c() #initialize a vector to contain opacity adjusted color to plot for allele frequency
		for (i in 1:nrow(g))
		{
		  tpc=c(tpc,opacitybyallelefreq(i,b,g))
		}

		#add row to table containing opacity adjusted color
		g$freqcol=tpc

		#plot an empty graph and add appropriately colored and sized line segments
		par(mai=c(2.25,1.6,0.8,0.4)) #origpar$mai=c(1.0,1.6,0.8,0.4)

		pdf(outfilename, width=7, height=4.35)

		plot(g$mhstart,g$pop,type="n",yaxt="n",xlim=c(min(g$mhstart),max(g$mhend)),
			ylim=c(min(g$pop)-1,max(g$pop)+1), xlab="", ylab="") #set up plot with no points showing
		axis(2, at=labelpos, labels=pools, las=1) 

		if (plottype==1)
		{
			#display segments with augmented ends, to fill empty space to next locus
			segments(g$mhstart,g$pop,g$mhnext,g$pop,col=g$freqcol,lwd=30,lend=1) #add colored lines of appropriate length for locus
		} else
		if (plottype==2)
		{
			#display microhaploblock segments with length corresponding to actual end points
			segments(g$mhstart,g$pop,g$mhend,g$pop,col=g$freqcol,lwd=30,lend=1) #add colored lines of appropriate length for locus
		}

		#annotate graph with functional regions
		if (gene=="BTC1L1")
		{
			fname="1" #name of feature
			f=c(17275,17590) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="2" #name of feature
			f=c(23104,23265) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="3" #name of feature
			f=c(23394,23527) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="4" #name of feature
			f=c(23783,23944) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="5" #name of feature
			f=c(24736,24916) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="6" #name of feature
			f=c(24999,25443) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="7" #name of feature
			f=c(25519,26335) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="8" #name of feature
			f=c(26761,27127) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])
		}
		dev.off()
	} #infilenames
} #genelist
	
par(mai=origpar$mai)
# End Main #

### END R MicrohaploblockHeatmap.r ###

### END PLOT MAJOR ALLELE MICROHAPLOBLOCK DIFFERENCES ###
### END BVBTC1 ALLELE MINING





### HS4 ALLELE MINING ###
### DETERMINE CONTIGS CONTAINING HS4 IN PATELLIFOLIA POOL ASSEMBLIES ###

module load blast+/2.9.0;

cd /home/pat.reeves/patellifolia/FlashedReadArchive;

#HS4
for i in $(seq 50 1 55);
  do blastn -db /home/pat.reeves/patellifolia/FlashedReadArchive/"$i"fraFinal/"$i"Hs1pro1REVdref.fasta -query /home/pat.reeves/patellifolia/seq/Hs4_Full_CDS.fa -out "$i"mapxHS4.txt;
  done;

#consolidated output from blast
#all searches above return x hit in procumbens/webbiana, y hits in all three patellaris because it is tetraploid

			Query= Rhomboid_Full_CDS

			Length=633
																				  Score        E
			pool    Sequences producing significant alignments:                          (Bits)     Value

			50      *jcf7180008065999                                                      425        4e-117
			50      *jcf7180008054617                                                      414        9e-114
			50      jcf7180008086815                                                      65.8       1e-08    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			50      jcf7180008063393                                                      65.8       1e-08    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			51      *jcf7180007596878                                                      425        4e-117
			51      *jcf7180007748729                                                      414        9e-114
			51      jcf7180007605177                                                      121        2e-25    only a 68 bp match, to query seq ATGGAGGCAGTATTCTGGATTATTCTTTTGAACTTTGTTATATACGGAGCAGAACACCTT
			51      jcf7180007746860                                                      65.8       1e-08    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			51      jcf7180007735918                                                      65.8       1e-08    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			52      *jcf7180008033456                                                      425        4e-117
			52      *jcf7180007998718                                                      409        4e-112
			52      *jcf7180008034409                                                      324        2e-86    matches two segments, ~350 bp combined, should be considered
			52      jcf7180008016940                                                      65.8       1e-08    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			52      jcf7180007949342                                                      65.8       1e-08    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			53      *jcf7180008836444                                                      425        3e-117
			53      jcf7180008863063                                                      65.8       8e-09    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			53      jcf7180008856006                                                      65.8       8e-09    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			54      *jcf7180009171055                                                      431        7e-119
			54      jcf7180009173230                                                      340        1e-91    matches 184bp exactly, contig is short, should be considered
			54      jcf7180009183662                                                      65.8       8e-09    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			54      jcf7180009162295                                                      65.8       8e-09    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			55      *jcf7180008572969                                                      431        6e-119
			55      jcf7180008540272                                                      65.8       7e-09    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG
			55      jcf7180008480775                                                      65.8       7e-09    only a 38 bp match, to query seq TACCAATTTGTGACCTCTGCATTTTGTCATTATAATTG

			        *use these for manual alignments and inspection using sequencher

#extract candidate contigs above from FRA for manual inspection
cd /home/pat.reeves/patellifolia/FlashedReadArchive;
sgrep '>jcf7180008065999_1_' 50fraFinal/50frasorted.fa | tr ' ' '\n' > 50jcf7180008065999.fa;
sgrep '>jcf7180008054617_1_' 50fraFinal/50frasorted.fa | tr ' ' '\n' > 50jcf7180008054617.fa;
sgrep '>jcf7180007596878_1_' 51fraFinal/51frasorted.fa | tr ' ' '\n' > 51jcf7180007596878.fa;
sgrep '>jcf7180007748729_1_' 51fraFinal/51frasorted.fa | tr ' ' '\n' > 51jcf7180007748729.fa;
sgrep '>jcf7180008033456_1_' 52fraFinal/52frasorted.fa | tr ' ' '\n' > 52jcf7180008033456.fa;
sgrep '>jcf7180007998718_1_' 52fraFinal/52frasorted.fa | tr ' ' '\n' > 52jcf7180007998718.fa;
sgrep '>jcf7180008034409_1_' 52fraFinal/52frasorted.fa | tr ' ' '\n' > 52jcf7180008034409.fa;
sgrep '>jcf7180008836444_1_' 53fraFinal/53frasorted.fa | tr ' ' '\n' > 53jcf7180008836444.fa;
sgrep '>jcf7180009171055_1_' 54fraFinal/54frasorted.fa | tr ' ' '\n' > 54jcf7180009171055.fa;
sgrep '>jcf7180008572969_1_' 55fraFinal/55frasorted.fa | tr ' ' '\n' > 55jcf7180008572969.fa;
mv 5[0-5]jcf*.fa /home/pat.reeves/patellifolia/AlleleMining/Hs4;
mv *mapxHS4.txt  /home/pat.reeves/patellifolia/AlleleMining/Hs4;

#Manual inspection of aligned sequences using sequencher showed shared indels can be used to
#distinguish orthologs.
#Within tetraploid patellaris one of the two sequences returned by blast is orthologous to
#the single sequence found in procumbens/webbiana



			*53jcf7180008836444  sense oriented but only 900 bp upstream
			*51jcf7180007596878  sense oriented but not full length over Hs4

			Hs4 locus 1 orthologs
			50jcf7180008065999:5622-9050
			51jcf7180007596878:1-3427
			52jcf7180008033456:5675-9103
			*53jcf7180008836444:1-3421
			54jcf7180009171055:28603-32037
			55jcf7180008572969:11767-15062
			* = proper sense orientation, use this one as single reference gene for mapping against

			Hs4 locus 2 orthologs
			50jcf7180008054617
			51jcf7180007748729
			52jcf7180007998718 + 52jcf7180008034409  these form a single (unassembled) contig spanning Hs4


#After identifying which contigs contain Hs4 locus 1 in each pool assembly, and the
#coordinates within those contigs that contain the Hs4 gene, extract the phased reads
#that map to that region from the FRA map. Just use samtools for this since the reads are
#already phased and quality filtered by virtue of being in the map
			
cd /home/pat.reeves/patellifolia/AlleleMining/Hs4;
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/50fraFinal/50fra.bam "jcf7180008065999:5622-9050" > 50Hs4L1.bam; #3429bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/51fraFinal/51fra.bam "jcf7180007596878:1-3427" > 51Hs4L1.bam; #3427bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/52fraFinal/52fra.bam "jcf7180008033456:5675-9103" > 52Hs4L1.bam; #3429bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/53fraFinal/53fra.bam "jcf7180008836444:1-3421" > 53Hs4L1.bam; #3421bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/54fraFinal/54fra.bam "jcf7180009171055:28603-32037" > 54Hs4L1.bam; #3435bp
samtools view -b -h /home/pat.reeves/patellifolia/FlashedReadArchive/55fraFinal/55fra.bam "jcf7180008572969:11767-15062" > 55Hs4L1.bam; #3296bp

#extract phased reads from bam files into a single fasta file
seq 50 1 55 | parallel --keep-order 'echo {}; samtools fasta {}Hs4L1.bam > {}Hs4L1.fa';

			50
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 699 reads
			51
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 297 reads
			52
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 585 reads
			53
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 787 reads
			54
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 813 reads
			55
			[M::bam2fq_mainloop] discarded 0 singletons
			[M::bam2fq_mainloop] processed 569 reads
			
#manual depth calculation from samtools stats 5[0-5]Hs4l1.finalaln.bam.TMP | grep 'sequences\|length' | grep -v ^#;
			pool	sequences	longest frag	total length	phased read count	phased read total length	cistron length	avg phased read depth
			50		699			24242			218121			698					193879						3429			56.54
			51		297			7781			89032			296					81251						3427			23.7
			52		585			23982			186469			584					162487						3429			47.39
			53		787			28536			242608			786					214072						3421			62.58
			54		813			49556			274714			812					225158						3435			65.55
			55		569			16110			158293			568					142183						3296			43.14

#map these phased reads to a single reference gene sequence (choose 53jcf7180008836444,
#which is in sense orientation and has 900 bp upstream of start codon)
#for Hs4L1 make 53jcf7180008836444_ref.txt from 53jcf7180008836444.fa, indexed with bwa index.
#you must rename the single reference sequence from "jcf7180008836444_1_28536" because a
#sequence by that name is included in the fasta files from above. It will cause hapx to
#fail in the calculation of $tmpf
cd /home/pat.reeves/patellifolia/AlleleMining/Hs4;
mkdir bam ref;
sed 's/jcf7180008836444_1_28536/53jcf7180008836444_ref/' 53jcf7180008836444.fa > ref/53jcf7180008836444_ref.txt;
bwa index ref/53jcf7180008836444_ref.txt;
pd=$(pwd);
t=40; #threads
for j in 50Hs4L1.fa 51Hs4L1.fa 52Hs4L1.fa 53Hs4L1.fa 54Hs4L1.fa 55Hs4L1.fa;
    do rg=$(cut -c1-2 <<<"$j"); #get the read group name 5[0-5]
      echo "$rg";
      readgroup="@RG\tID:$rg\tSM:$rg\tPL:illumina\tLB:na\tPU:na";
      bwa mem -R "$readgroup" -t "$t" "$pd"/ref/53jcf7180008836444_ref.txt "$pd"/"$j" | \
                  samtools sort -O BAM --threads "$t" | \
                  samtools view -F 2048 -O BAM > "$pd"/bam/"$rg"Hs4l1.finalaln.bam.TMP; #-F 2048 excludes supplementary alignments
                  samtools index "$pd"/bam/"$rg"Hs4l1.finalaln.bam.TMP;
    done;
sambamba merge -t40 -p "$pd"/bam/Hs4l1.finalaln.merge.bam "$pd"/bam/5[0-5]Hs4l1.finalaln.bam.TMP;


### END DETERMINE CONTIGS CONTAINING HS4 IN PATELLIFOLIA POOL ASSEMBLIES ###


### OPTIONAL: COUNT TOTAL, DISTINCT, AND PRIVATE ALLELES ###
#rsync results from HPC ceres to HPC blip so you can get access to multi-node gnu parallel

#run a final hapx to count microhaplotypes at each position (takes ~25 min)
cd /share/space/reevesp/patellifolia/AlleleMining/Hs4; #go here just to get the working directory specified
pd=$(pwd);
time hapx.sh -r "$pd"/ref/53jcf7180008836444_ref.txt \
        -b "$pd"/bam/Hs4l1.finalaln.merge.bam \
        -o Hs4l1Finalhapx \
        -f 0 -q 1 -x -ssh /home/reevesp/machines \
        -s <(for i in $(seq 1 1 4000); do echo 53jcf7180008836444_ref:"$i"-"$i"; done;); #expand window around locus in ref somewhat (1-3421)

#postprocess haploblock counts in log.txt into something plottable
cd /share/space/reevesp/patellifolia/AlleleMining/Hs4/Hs4l1Finalhapx;

### BEGIN BASH ###
mytd1() {
        i=$1; #a position on the reference is incoming
        a=$(grep "$i"\\."$k" "$pd"/numblocks.txt | cut -d: -f4); #number of unique haploblock read pairs mapped to position i in the readgroup
        if [[ "$a" == "" ]]; then a="?"; fi; #if no data at position i set number of haploblock read pairs to ?
        echo "$i"."$k $a"; #report result to parallel statement
}
export -f mytd1;

mytd2() {
        i=$1; #a position on the reference is incoming
        b=$(grep "$i"\\."$k" "$pd"/numblocks.txt | cut -d$'\t' -f2 | cut -d: -f1); #total number haploblock read pairs mapped to position i in the readgroup
        if [[ "$b" == "" ]]; then b="?"; fi; #if no data at position i set to ?
        echo "$i"."$k $b"; #report result to parallel statement
}
export -f mytd2;

mypa() {
       i=$1; #contig:site-range
       j=$(grep "$i" "$pd"/counts.txt | grep -v global);
       names=$(cut -d$'\t' -f1 <<<"$j" | sed 's/$/ 0/');
       counts=$(cut -d$'\t' -f2 <<<"$j");
       nc=$(head -1 <<<"$counts" | awk -F: '{print NF}'); #number of alleles
       nr=$(wc -l <<<"$counts"); #number of read groups
       nz=$(( $nr - 1 )); #number of zeroes needed in a column of allele counts for a private allele to exist
       for m in $(seq 1 1 $nc);
         do col=$(cut -d: -f$m <<<"$counts"); #extract the column of allele counts
           fl=$(grep -n -v 0 <<<"$col"); #show line numbers where alleles are present
           
           #if only one line of the column of data has alleles, it is a private allele
           if [[ $(echo "$fl" | wc -l) == 1 ]];
           then r=$(cut -d: -f1 <<<"$fl"); #row number of private allele
             names=$(awk -F' ' -v r=$r '{if (NR==r) $2++; print}' <<<"$names"); #index up by one the row with the private allele
           fi;
         done;
         
       #report result to parallel statement
       echo "$names";
}
export -f mypa;

mypadpa() {
          aa=$1;
          cc=$(grep "$aa.$dd" "$pd"/names.txt);
          if [[ "$cc" == "" ]];
          then echo "$aa.$dd ?"; #if readgroup not found print a ?
          else echo "$cc";
          fi;
}
export -f mypadpa;


pd=$(pwd); export pd;

#count total and distinct haploblocks
grep ^# log.txt | tr '-' '_' > numblocks.txt;
a=$(rev numblocks.txt | cut -d. -f2- | rev | uniq); #list of contig:site-ranges
b=$(rev numblocks.txt | cut -d. -f1 | rev | cut -d$'\t' -f1 | sort -u); export b; #determine set of possible read groups

#iterate over read groups to count haploblocks
for j in $b;
  do k="$j"; export k; #do this so iterator can be sent to nodes using --sshloginfile
    echo "$j: counting distinct alleles";
    >"$pd"/"$j".distallel.txt; #initialize output file with a header for count of distinct alleles within populations
    echo "$a" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env mytd1 --env pd --env k /home/reevesp/bin/parallel --jobs 96 --env mytd1 --env pd --env k mytd1 >> "$pd"/"$j".distallel.txt;
    #echo "$a" | parallel --jobs 1 --pipe -N960 --env mytd1 --env pd --env j /home/reevesp/bin/parallel --jobs 96 --env mytd1 --env pd --env j mytd1 >> "$pd"/"$j".distallel.txt;
    
    echo "$j: counting all alleles";
    >"$pd"/"$j".totallel.txt; #initialize output file with a header for count of all alleles within populations
    echo "$a" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env mytd2 --env pd --env k /home/reevesp/bin/parallel --jobs 96 --env mytd2 --env pd --env k mytd2 >> "$pd"/"$j".totallel.txt;
  done;


#count private alleles
grep ^'@' log.txt | tr '-' '_' | cut -d$'\t' -f1,3 | sort -t_ -k2,2n > freqs.txt;
grep ^'@' log.txt | tr '-' '_' | cut -d$'\t' -f1,2 | sort -t_ -k2,2n > counts.txt;

c=$(rev counts.txt | cut -d. -f2- | rev | uniq); #list of contig:site-ranges
>names.txt;
echo "$c" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env pd --env mypa /home/reevesp/bin/parallel --jobs 96 --env pd --env mypa mypa >> names.txt; #counting sub

#pad the file containing counts of private alleles with "?" for contig:site-range_readgroup combinations that weren't found
for bb in $b;
do >"$bb".privallel.txt; #create an output file for the readgroup
  dd="$bb"; export dd; #transfer iterator to variable that can be exported for gnu parallel
  echo "$a" | sed 's/#/@/g' | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env pd --env dd --env mypadpa \
                              /home/reevesp/bin/parallel --jobs 96 --env pd --env dd --env mypadpa mypadpa >> "$bb".privallel.txt;
done;

#sort and tab delimit output
for i in $b;
  do sort -t'_' -k1,1 "$i".distallel.txt | sort -t'_' -k2,2n | sed 's/^#//' | sed "1i position $i.distallel" | tr ' ' '\t' > "$i".2.distallel.txt;
    sort -t'_' -k1,1 "$i".totallel.txt | sort -t'_' -k2,2n | sed 's/^#//' | sed "1i position $i.totallel" | tr ' ' '\t'  > "$i".2.totallel.txt;
    sort -t'_' -k1,1 "$i".privallel.txt | sort -t'_' -k2,2n | sed 's/^@//' | sed "1i position $i.privallel" | tr ' ' '\t'  > "$i".2.privallel.txt;
  done;
#swap back to original filename
for i in $b; 
  do mv "$i".2.distallel.txt "$i".distallel.txt;
    mv "$i".2.totallel.txt "$i".totallel.txt;
    mv "$i".2.privallel.txt "$i".privallel.txt;
  done;
### END BASH ###


#consolidate files into 1, this part depends on number of read groups, is not universal
#you can tailor the columns cut depending on the number of read groups included
#for <= 6 read groups, it turns out the same
paste -d$'\t' global.distallel.txt RG_Z_*.distallel.txt global.totallel.txt RG_Z_*.totallel.txt global.privallel.txt RG_Z_*.privallel.txt \
  | cut -d$'\t' -f1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42 | sed 's/\.global//'> summary.txt;

#Pools differed in total number of reads. Scale proportionally by total read count relative to pool 51,
#which had the fewest
			L1
			pool:cols:scalar
			50:3,10,17:1.71
			51:4,11,18:1.00
			52:5,12,19:1.23
			53:6,13,20:1.83
			54:7,14,21:1.18
			55:8,15,22:1.20
			L1global
			2,9,16
			$3..$8=$2
			$10..$15=$9
			$17..$22=$16

cd /share/space/reevesp/patellifolia/AlleleMining/Hs4/Hs4l1Finalhapx;

summ=$(tail -n +2 summary.txt | tr "\t" " "); #acquire summary.txt as a variable
h=$(head -1 summary.txt);
#adjust values for each read group using scaling values above
l1="50:3,10,17:1.71 51:4,11,18:1.00 52:5,12,19:1.23 53:6,13,20:1.83 54:7,14,21:1.18 55:8,15,22:1.20";
for i in $l1;
  do echo "$i";
    c=$(echo "$i" | cut -d: -f2); #columns to scale
    s=$(echo "$i" | cut -d: -f3); #scalar, divide by this
    #cycle through list of columns to scale
    for cc in $(echo "$c" | sed 's/,/ /g');
      do summ=$(echo "$summ" | awk -F$' ' -v cc=$cc -v s=$s '{$cc=$cc/s; print $0}');
      done;
  done;

#adjust the so-called 'global' value, which is just the sum of the read group values
l1global="2:3:8 9:10:15 16:17:22"; #columns to sum
for i in $l1global;
  do s=$(echo "$i" | cut -d: -f1); #column to put sum in
    st=$(echo "$i" | cut -d: -f2); #first column to sum
    se=$(echo "$i" | cut -d: -f3); #last column to sum
    summ=$(echo "$summ" | awk -F' ' -v s=$s -v st=$st -v se=$se '{for(j=st;j<=se;j++)x+=$j;$s=x; print $0; x=0}');
  done;

#write
echo "$h" > summaryscaled.txt;
echo "$summ" | awk -F' ' '{for(i=1; i<=NF; i++) {if($i=="0") $i="?"}; print $0}' | tr ' ' '\t' >> summaryscaled.txt;

### END OPTIONAL: COUNT TOTAL, DISTINCT, AND PRIVATE ALLELES ###



### IDENTIFY MICROHAPLOBLOCKS ACROSS HS4 ###

#Use hapxm.sh to find microhaploblocks in the bam files for each EL10 sco output by hapx -mb.
#By using -db and -va you can produce an
#optional file, mhendsvar.txt, that contains a variant optimized tiling path across the locus thru the
#microhaploblocks, which can be used directly as input into -u, to specify the universal microhaploblock set.
#hapxm option -va calculates a tiling path that collects the most variable (then, for ties, 
#most depth, then longest) microhaploblocks at each site.

#Hs4l1, ~4 minutes on all nodes
#determine range to consider
cd /share/space/reevesp/patellifolia/AlleleMining/Hs4/bam;
samtools view Hs4l1.finalaln.merge.bam | grep -v ^'jcf' | awk -F$'\t' '$4!=0{print $0}' | cut -d$'\t' -f4 | sort -n | head; #leftmost mapping position = 1
samtools view Hs4l1.finalaln.merge.bam | grep -v ^'jcf' | awk -F$'\t' '$4!=0{print $0}' | cut -d$'\t' -f4 | sort -n | tail; #rightmost mapping position = 4830 + 150 = 4980

cd /share/space/reevesp/patellifolia/AlleleMining/Hs4;
mkdir hapxmHs4L1;
cd hapxmHs4L1;
time hapxm.sh -b /share/space/reevesp/patellifolia/AlleleMining/Hs4/bam/Hs4l1.finalaln.merge.bam \
            -ssh /home/reevesp/machines \
            -o hxm1 -db -va -s <(for i in $(seq 1 1 4980); do echo 53jcf7180008836444_ref:"$i"-"$i"; done;)

#use hapxm.sh to find microhaploblocks in the finalaln.bam.TMP file from each read group using the universal
#microhaplotypes discovered above. ~12 minutes on all nodes
time for i in $(seq 50 1 55);
  do echo $i;
    hapxm.sh -b /share/space/reevesp/patellifolia/AlleleMining/Hs4/bam/"$i"Hs4l1.finalaln.bam.TMP \
            -ssh /home/reevesp/machines \
            -o hxm"$i"onhxm1 -u /share/space/reevesp/patellifolia/AlleleMining/Hs4/hapxmHs4L1/hxm1/mhendsvar.txt \
            -s <(for i in $(seq 1 1 4980); do echo 53jcf7180008836444_ref:"$i"-"$i"; done;)
  done;
### END IDENTIFY MICROHAPLOBLOCKS ACROSS HS4 ###



### FIND MICROHAPLOBLOCKS WHERE MAJOR ALLELE DIFFERS BETWEEN POOLS ###

cd /share/space/reevesp/patellifolia/AlleleMining/Hs4/hapxmHs4L1;
rhead=$(grep ^"#" hxm50onhxm1/hapxmlog.txt | cut -d' ' -f1-4);
rhead=$(tail -n +2 <<<"$rhead");#remove col head from row heads

outmajall="";
outmajfreq="";
for i in $(seq 50 1 55);
  do echo $i;
    majall=$(grep ^"#" hxm"$i"onhxm1/hapxmlog.txt | tail -n +2 | cut -d' ' -f9 | cut -d: -f1); #extract the allele with highest frequency
    majfreq=$(grep ^"#" hxm"$i"onhxm1/hapxmlog.txt | tail -n +2 | cut -d' ' -f8 | cut -d: -f1); #extract the frequency of the major allele
    outmajall=$(paste -d: <(echo "$outmajall") <(echo "$majall"));
    outmajfreq=$(paste -d: <(echo "$outmajfreq") <(echo "$majfreq"));
  done;
outmajall=$(sed 's/^://' <<<"$outmajall"); #remove leading colon
outmajfreq=$(sed 's/^://' <<<"$outmajfreq"); #remove leading colon

#recode mhblocks from DNA sequence to integers
set -f; #you have to undo globbing before running this command since some values are solely asterisks
outrecode="";
while read l;
do m=$(echo "$l" | tr ':' '\n' | sed 's/^$/0/' | tr '\n' ':' | sed 's/:$//'); #substitute 0 when no allele was found
  a=$(echo "$m" | tr ':' '\n' | sort -u | grep -v 0); #list containing DNA sequences of alleles, excluding 0 (missing)
  j=1; #numeric allele name
  ll="$m";
  #if [[ "$a" == "" ]];
  #then outrecode+="$l"$'\n'; #deal with mhloci with no alleles by placing :::::
  #else
    for aa in $a; 
      do ll=${ll//"$aa"/"$j"}; #built-in bash replace obviates need to escape asterisks using sed
        j=$(( $j + 1 ));
      done;
     outrecode+="$ll"$'\n';
  #fi;
done <<<"$outmajall";
outrecode=$(sed '/^$/d' <<<"$outrecode"); #remove terminal blank line
set +f; #redo globbing

#assemble output
outp1=$(paste -d' ' <(echo "$rhead") <(echo "$outmajall") <(echo "$outrecode") <(echo "$outmajfreq")); #add row header to major allele list

#find all mh loci with different major alleles, label those rows with @
outp2="";
outp2=$(while read l;
do a=$(echo "$l" | cut -d' ' -f5 | tr ':' '\n'| sed '/^$/d' | sort | uniq | wc -l);
  if [[ "$a" > 1 ]];
  then echo "@$l";
  else echo "$l"
  fi;
done <<<"$outp1";
);

#write out the summary
echo "$outp2" > hxmsummary.txt;

### END FIND MICROHAPLOBLOCKS WHERE MAJOR ALLELE DIFFERS BETWEEN POOLS ###



### CREATE INPUT FILES FOR R ROUTINE MicrohaploblockHeatmap.r ###

cd /share/space/reevesp/patellifolia/AlleleMining/Hs4/hapxmHs4L1;
c=$(cut -d' ' -f2-8 hxmsummary.txt); #get mhstart mhend mhlength seqs majalleles freqs for all mh majalleles
d=$(grep ^@ hxmsummary.txt | cut -d' ' -f2-8); #get mhstart mhend mhlength seqs majalleles freqs for mh majalleles that differ between pops

for k in "$c" "$d";
  do hdr="pop mhstart mhend length majallele freq seq";
    e="";
    while read l;
      do sqc=$(echo "$l" | cut -d' ' -f4); #allele sequence column
        mac=$(echo "$l" | cut -d' ' -f5); #major alleles column
        frc=$(echo "$l" | cut -d' ' -f6); #frequencies column
        j=1; #j indexes position in horizontal allele calls and freqs
        for i in $(seq 50 1 55);
          do rhead=$(echo -n "$i ";echo "$l" | cut -d' ' -f1-3); #row header labeled by population name
            sq=$(echo "$sqc" | cut -d: -f$j); #allele sequence for population $i
            ma=$(echo "$mac" | cut -d: -f$j); #major allele for population $i
            fr=$(echo "$frc" | cut -d: -f$j); #major allele frequency for population $i
            e+="$rhead $ma $fr $sq"$'\n';
    
            j=$(( $j + 1 ));
          done;
      done<<<"$k";
    e=$(sed '/^$/d' <<<"$e");
    e="$hdr"$'\n'"$e";
    
    
    if [[ "$k" == "$c" ]];
    then echo "$e" | tr ' ' '\t' > hxmsummaryRinAll.txt;
    elif [[ "$k" == "$d" ]]; 
    then echo "$e" | tr ' ' '\t' > hxmsummaryRinDiff.txt;
    fi;
  done;

### END CREATE INPUT FILES FOR R ROUTINE MicrohaploblockHeatmap.r ###



### CALCULATE DESCRIPTIVE STATISTICS FOR MICROHAPLOBLOCKS IN VARIANT-RICH TILING ARRAY ###

#myst() calculates basic stats
myst() {
  awk '{
    sum = 0;
    M = 0;
    S = 0;
    for (k=1; k <= NF; k++) {
      sum += $k;
      x = $k;
      oldM = M;
      M = M + ((x - M)/k);
      S = S + (x - M)*(x - oldM);
    }
    var = S/(NF - 1);
    print "n=" NF " mean=" sum/(NF) " var=" var " sd=" sqrt(var);
  }' <<< $1
}
export myst;

cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/Hs4/hapxmHs4L1;

#Choose boundaries by visual inspection of hxmsummaryRinAll.txt and hxmsummaryRinDiff.txt in order to
#trim off long mhs and poor alignment regions from ends, which are artifacts. Choice is subjective.

for i in {50..55};
  do a=$(sed -n '/^#53jcf7180008836444_ref 491 492/,$p' hxm"$i"onhxm1/hapxmlog.txt | \
         sed -e '/^#53jcf7180008836444_ref 3740 3740/,$d'); #sed -n $p prints lines after including match, -e $d deletes lines after including match
     b=$(echo "$a" | awk -F' ' '$4!=1{print $0}'); #exclude loci with mhlength=1
     #b=$(echo "$b" | awk -F' ' '$6!=1{print $0}'); #exclude loci with only 1 variant (don't use this)
     nmh=$(echo "$b" | wc -l); #number of mhs across minimum tiling path
     
     #calc stats
     mhl=$(echo "$b" | cut -d' ' -f4 | tr "\n" " " | sed 's/ $//'); #stats on mhlength
     mhlengthst=$(myst "$mhl");
     ns=$(echo "$b" | cut -d' ' -f5 | tr "\n" " " | sed 's/ $//'); #stats on numseq
     numseqst=$(myst "$ns");
     na=$(echo "$b" | cut -d' ' -f6 | tr "\n" " " | sed 's/ $//'); #stats on numalleles
     numallelesst=$(myst "$na");
     
     #calc mean length of mhs
     #c=$(echo "$b" | awk -F' ' '{s+=$4} END {print s}'); #sum of mhlength
     #ml=$(echo "scale=4; $c / $nmh" | bc);
     
     echo "$i $mhlengthst | $numseqst | $numallelesst"
  done;

#result, basic stats, mh with 1 variant included:
			50 n=754 mean=3.79045 var=6.64926 sd=2.57862 | n=754 mean=30.7692 var=74.9852 sd=8.6594 | n=754 mean=1.9695 var=0.435987 sd=0.660293
			51 n=754 mean=3.79045 var=6.64926 sd=2.57862 | n=754 mean=9.67905 var=58.5688 sd=7.65303 | n=754 mean=1.01326 var=0.719611 sd=0.848299
			52 n=754 mean=3.79045 var=6.64926 sd=2.57862 | n=754 mean=28.6645 var=42.4118 sd=6.51244 | n=754 mean=1.94032 var=0.393512 sd=0.627305
			53 n=754 mean=3.79045 var=6.64926 sd=2.57862 | n=754 mean=38.2294 var=71.9274 sd=8.481 | n=754 mean=2.47082 var=0.79928 sd=0.894025
			54 n=754 mean=3.79045 var=6.64926 sd=2.57862 | n=754 mean=35.3064 var=51.0707 sd=7.14638 | n=754 mean=2.23475 var=0.557038 sd=0.74635
			55 n=754 mean=3.79045 var=6.64926 sd=2.57862 | n=754 mean=26.565 var=43.0057 sd=6.55788 | n=754 mean=1.62334 var=0.447582 sd=0.669016


#calculate the proportion of loci containing only indel variants
for i in {50..55};
  do 
     a=$(sed -n '/^#53jcf7180008836444_ref 491 492/,$p' hxm"$i"onhxm1/hapxmlog.txt | \
         sed -e '/^#53jcf7180008836444_ref 3740 3740/,$d'); #sed -n $p prints lines after including match, -e $d deletes lines after including match
     b=$(echo "$a" | awk -F' ' '$4!=1{print $0}'); #exclude lines with mhlength=1
     c=$(echo "$b" | awk -F' ' '$6!=1{print $0}' | cut -d' ' -f9); #collect variants, excluding loci with only 1 variant
     #c=$(echo "$b" | cut -d' ' -f9); #collect variants, including loci with only 1 variant
     
     #count loci that only contain indel variants
     indelloc=0; #number of loci that contain only indel type variants
     snploc=0
     while read d;
       do e=$(echo "$d" | tr ':' '\n' | awk NF=NF FS= | rs -c' ' -C' ' -T | sed 's/ //g' \
              | sed 's/\*//g' | grep -v ^"A\+"$ | grep -v ^"C\+"$ | grep -v ^"G\+"$ | grep -v ^"T\+"$ \
              | wc -l); #split to rows, add space btw each character, transpose, remove space delimiter
                        #delete asterisks, remove lines that are all repeated ACG or T
                        #count number of lines remaining (if >0 the locus contains at least 1 SNP variant
                        #if =0 the locus contains only indel variants)
          if (( $e == 0 )); then (( indelloc++ )); fi;
          if (( $e > 0 )); then (( snploc++ )); fi;
       done <<<"$c";
       
       totloc=$(( indelloc + snploc ));
       indelpro=$(echo "scale=4; $indelloc / $totloc" | bc);
       snppro=$(echo "scale=4; $snploc / $totloc" | bc);
     
     echo "$i $indelloc $indelpro $snploc $snppro $totloc";
  done;

#results, single variant loci excluded:
			50 591 .9784 13 .0215 604
			51 469 .9894 5 .0105 474
			52 577 .9763 14 .0236 591
			53 600 .8836 79 .1163 679
			54 609 .9076 62 .0923 671
			55 367 .9314 27 .0685 394

#repeat calculating the proportion of loci containing only indel variants for the combined dataset
cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/Hs4/hapxmHs4L1;
for i in hxm1;
  do 
     a=$(sed -n '/^#53jcf7180008836444_ref 491 492/,$p' "$i"/hapxmlogvar.txt | \
         sed -e '/^#53jcf7180008836444_ref 3740 3740/,$d'); #sed -n $p prints lines after including match, -e $d deletes lines after including match
     b=$(echo "$a" | awk -F' ' '$4!=1{print $0}'); #exclude lines with mhlength=1
     c=$(echo "$b" | awk -F' ' '$6!=1{print $0}' | cut -d' ' -f9); #collect variants, excluding loci with only 1 variant
     #c=$(echo "$b" | cut -d' ' -f9); #collect variants, including loci with only 1 variant
     
     #count loci that only contain indel variants
     indelloc=0; #number of loci that contain only indel type variants
     snploc=0
     while read d;
       do e=$(echo "$d" | tr ':' '\n' | awk NF=NF FS= | rs -c' ' -C' ' -T | sed 's/ //g' \
              | sed 's/\*//g' | grep -v ^"A\+"$ | grep -v ^"C\+"$ | grep -v ^"G\+"$ | grep -v ^"T\+"$ \
              | wc -l); #split to rows, add space btw each character, transpose, remove space delimiter
                        #delete asterisks, remove lines that are all repeated ACG or T
                        #count number of lines remaining (if >0 the locus contains at least 1 SNP variant
                        #if =0 the locus contains only indel variants)
          if (( $e == 0 )); then (( indelloc++ )); fi;
          if (( $e > 0 )); then (( snploc++ )); fi;
       done <<<"$c";
       
       totloc=$(( indelloc + snploc ));
       indelpro=$(echo "scale=4; $indelloc / $totloc" | bc);
       snppro=$(echo "scale=4; $snploc / $totloc" | bc);
     
     echo "$i $indelloc $indelpro $snploc $snppro $totloc";
  done;


#results, single variant loci excluded:
			hxm1 605 .8045 147 .1954 752

#calculate major variant frequency
cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/Hs4/hapxmHs4L1;
a=$(sed -n '/^#53jcf7180008836444_ref 491 492/,$p' hxmsummary.txt | \
    sed -e '/^#53jcf7180008836444_ref 3740 3740/,$d'); #sed -n $p prints lines after including match, -e $d deletes lines after including match
b=$(echo "$a" | awk -F' ' '$4!=1{print $0}'); #exclude lines with mhlength=1
c=$(echo "$b" | cut -d' ' -f7); #collect major variant frequencies, including loci with only 1 variant

#calc stats
for i in {50..55};
  do j=$(( $i - 49 )); #calculate a column index for each pool
    vf=$(echo "$c" | cut -d: -f"$j" | grep -v ^0$ | tr "\n" " " | sed 's/ $//'); #get variant freqs, remove zeroes (no variant called)
    varfreqs=$(myst "$vf");
    echo "$i $varfreqs";
  done;
			50 n=754 mean=0.898217 var=0.00558517 sd=0.074734
			51 n=511 mean=0.924333 var=0.00883043 sd=0.0939704
			52 n=754 mean=0.895995 var=0.00610176 sd=0.0781138
			53 n=753 mean=0.863518 var=0.0116013 sd=0.10771
			54 n=754 mean=0.868737 var=0.00964737 sd=0.098221
			55 n=754 mean=0.937211 var=0.0058436 sd=0.0764434

### END CALCULATE DESCRIPTIVE STATISTICS FOR MICROHAPLOBLOCKS IN VARIANT-RICH TILING ARRAY ###



### PLOT MAJOR ALLELE MICROHAPLOBLOCK DIFFERENCES ###

#rsync everything back to local machine for work in R
cd /Users/wichita/Desktop/telework/patellifolia/AlleleMining/Hs4;

### BEGIN R MicrohaploblockHeatmap.r ###
#install.packages("RColorBrewer")

options(error = recover)
rm(list=ls()) 
library("RColorBrewer")

# Functions #
	opacitybyallelefreq <- function(i,b,g)
	{
	  rr=rgb(b[1,i],b[2,i],b[3,i],alpha=g$freq[i]*255,max=255)
	  return(rr)
	}

	plotfeatures=function(f,fname,bbegin,eend)
	{
	  segments(f[1],eend+1,f[2],eend+1,col="black",lwd=20,lend=1)
	  segments(f[1],bbegin-1,f[1],eend+1,col="black",lwd=1,lend=1)
	  segments(f[2],bbegin-1,f[2],eend+1,col="black",lwd=1,lend=1)
	  text(f[1]+((f[2]-f[1])/2),eend+1,fname,col="white",cex=0.5,font=1)
	}
# End Functions #

# Main #
origpar=par() #gather initial settings
plottype=1 #1, fill empty space to next mh locus; 2, plot actual end points of mh

genelist=c("Hs4L1")
for (gene in genelist)
{
	print(gene)
	
	if (gene=="Hs4L1")
	{
		setwd("/Users/wichita/Desktop/telework/patellifolia/AlleleMining/Hs4/hapxmHs4L1")
		pools=c("pat S", "pat T", "pat S", "pro H", "pro T", "web GC") #labels in order 50:55 (Hs4L1)
		labelpos=50:55
	} 	

	#import data table
	infilenames=c("hxmsummaryRinAll.txt","hxmsummaryRinDiff.txt")
	for (infilename in infilenames)
	{
		print(infilename)
		
		outfilename=paste(gene,gsub("hxmsummaryRin","",infilename),sep="")
		outfilename=gsub(".txt",".pdf",outfilename)

		g <- read.table(infilename, header=TRUE, sep="\t")

		#add a column that defines the start of the next microhaploblock allele for the current one
		nd=length(unique(g$pop)) #number of lines to delete from mhstart to shift it for mhnext
		shft=g$mhstart[(nd+1):length(g$mhstart)] #remove first nd elements from mhstart
		shft=c(shft,tail(g$mhend,nd))
		g$mhnext=shft #add column that will extend the stop point to the beginning of the next locus, for segment length

		#create a list of hex rgb colors with opacity proportional to allele frequency
		#rcbcolors=col2rgb(brewer.pal(n = 8, name = "Dark2")) #get hex codes for 8 RColorBrewer colors, convert to rgb
		#rcbcolors=cbind(col2rgb("black"), rcbcolors) #add black in position 1
		#b=rcbcolors[,g$majallele+1] #get the RColorBrewer Dark2 color for all major alleles, as defined by integer
		b=col2rgb(g$majallele+1) #color by major allele state, add 1 so that missing data (0) means black (r color value 1)
								 #this only works until there are more than 7 colors (r has only 8 numerically coded numbers)
		tpc=c() #initialize a vector to contain opacity adjusted color to plot for allele frequency
		for (i in 1:nrow(g))
		{
		  tpc=c(tpc,opacitybyallelefreq(i,b,g))
		}

		#add row to table containing opacity adjusted color
		g$freqcol=tpc

		#plot an empty graph and add appropriately colored and sized line segments
		par(mai=c(2.25,1.6,0.8,0.4)) #origpar$mai=c(1.0,1.6,0.8,0.4)

		pdf(outfilename, width=7, height=4.35)

		plot(g$mhstart,g$pop,type="n",yaxt="n",xlim=c(min(g$mhstart),max(g$mhend)),
			ylim=c(min(g$pop)-1,max(g$pop)+1), xlab="", ylab="") #set up plot with no points showing
		axis(2, at=labelpos, labels=pools, las=1) 

		if (plottype==1)
		{
			#display segments with augmented ends, to fill empty space to next locus
			segments(g$mhstart,g$pop,g$mhnext,g$pop,col=g$freqcol,lwd=30,lend=1) #add colored lines of appropriate length for locus
		} else
		if (plottype==2)
		{
			#display microhaploblock segments with length corresponding to actual end points
			segments(g$mhstart,g$pop,g$mhend,g$pop,col=g$freqcol,lwd=30,lend=1) #add colored lines of appropriate length for locus
		}

		#annotate graph with functional regions
		if (gene=="Hs4L1")
		{
			fname="1" #name of feature
			f=c(897,962) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="2" #name of feature
			f=c(1633,1718) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="3" #name of feature
			f=c(1835,1875) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="4" #name of feature
			f=c(1984,2240) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])

			fname="5" #name of feature
			f=c(2739,2921) #position of feature
			plotfeatures(f,fname,labelpos[1],labelpos[length(labelpos)])
		}
		dev.off()
	} #infilenames
} #genelist
	
par(mai=origpar$mai)
# End Main #

### END R MicrohaploblockHeatmap.r ###

### END PLOT MAJOR ALLELE MICROHAPLOBLOCK DIFFERENCES ###
### END HS4 ALLELE MINING













