# sam2ngs
Process the sequence after the genome and reads are aligned.
## Installation
<pre><code>
wget -c https://github.com/zxgsy520/sam2ngs/archive/v1.1.0.tar.gz
tar -zxvf v1.1.0.tar.gz
cd sam2ngs-1.1.0
 chmod 755 *
</code></pre>
or
<pre><code>
git clone https://github.com/zxgsy520/sam2ngs.git
cd sam2ngs
chmod 755 *
</code></pre>
## Instructions
sam2ngs: Extract the sequence on the genome and the second-generation data alignment (sam).
<pre><code>
usage: sam2ngs [-h] [-i FILE] -g FILE [-d STR] [-p STR]

name:
    sam2ngs Select no reads on the map

attention:
    sam2ngs -i input.sam -g genome.fa -p name
    minimap2 -t 8 -ax sr genome.fa r1.fq r2.fq|samblaster -a |sam2ngs -g genome.fa -p name

version: 1.1.0
contact:  Xingguo Zhang <113178210@qq.com>        

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        Input files in sam.
  -g FILE, --genome FILE
                        Inport genome files.
  -d STR, --depth STR   Set the reads depth selected by each contig(Select all=all), default=200
  -p STR, --prefix STR  Inport the file prefix.
  </code></pre>
  filter_paf：Filter low-quality comparison data (paf) of genome and reads.
  <pre><code>
  usage: filter_paf [-h] [-i FILE] [-id FLOAT] [-mq INT]

name:
    filter_paf.py  Filter the reads comparison result file

attention:
    filter_paf.py -i file.paf >file_clean.paf
    cat file.paf | filter_paf.py >file_clean.paf

version: 1.1.0
contact:  Xingguo Zhang <113178210@qq.com>        

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        Input files in paf
  -id FLOAT, --identities FLOAT
                        Set the minimum identities of reads comparison, default=50.0
  -mq INT, --mapq INT   Set the minimum quality value for reads alignment, default=2
   </code></pre>
   sam2unmap：Extract the unaligned sequence of the genome and reads．
   <pre><code>
   usage: sam2unmap [-h] [-i FILE] [-mq INT]

name:
    sam2unmap Recover unmap sequence

attention:
    sam2unmap -i input.sam >unmap.fastq
    minimap2 -t 8 -ax sr genome.fa r1.fq r2.fq|sam2unmap >unmap.fastq

version: 1.1.0
contact:  Xingguo Zhang <113178210@qq.com>        

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        Input files in sam.
  -mq INT, --minmapq INT
                        Set the minimum quality value for reads alignment, default=20
 </code></pre>
   
