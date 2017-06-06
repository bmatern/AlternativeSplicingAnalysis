# This file is part of AlternativeSplicingAnalysis.
#
# AlternativeSplicingAnalysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# AlternativeSplicingAnalysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with AlternativeSplicingAnalysis. If not, see <http://www.gnu.org/licenses/>.



source /home/ben/anaconda2/bin/activate minionvironment

cd /home/ben/Github/AlternativeSplicingAnalysis/src

#python AnalyzeAlternativeSplicing.py \
#   --reads=/home/minion/MinIONData/cDNA_C04_GMAP/22677cDNA.fastq \
#   --reference=/home/minion/MinIONData/cDNA_C04_GMAP/C04010101FullSequence.fasta \
#    --odir=/home/minion/MinIONData/cDNA_C04_GMAP/splicingresults \
#    --exonmap=/home/minion/MinIONData/cDNA_C04_GMAP/C04010101.Features.txt

python AnalyzeAlternativeSplicing.py \
    --reads=/home/minion/MinIONData/cDNAReads_Oct4.2016_TimoExperiment/extracts/cDNAOct4_2017_unbarcoded_2D_reads.fastq \
    --reference=/home/minion/MinIONData/cDNAReads_Oct4.2016_TimoExperiment/inputfiles/C04010101FullSequence.fasta \
    --odir=/home/minion/MinIONData/cDNAReads_Oct4.2016_TimoExperiment/splicingresults \
    --exonmap=/home/minion/MinIONData/cDNAReads_Oct4.2016_TimoExperiment/inputfiles/C04010101.Features.txt
    
#--reads=/home/minion/MinIONData/cDNAReads_Oct4.2016_TimoExperiment/extracts/cDNAOct4_2017_unbarcoded_2D_reads.fastq \
#--reads=/home/minion/MinIONData/cDNAReads_Oct4.2016_TimoExperiment/inputfiles/22677cDNA.Short.fastq \

source /home/ben/anaconda2/bin/deactivate
