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


# Set some variables, like the pros do.


#9.  SAM output format
#=====================

#GSNAP can generate SAM output format, by providing the "-A sam" flag
#to GSNAP.  In addition, GMAP can also print its alignments in SAM
#output, using the "-f samse" or "-f sampe" options, for single-end or
#paired-end data.  The sampe option will generate SAM flags to indicate

ExonsOnlyReferenceInputFile="C04010101ExonsOnly.fasta"
FullGeneReferenceInputFile="C04010101FullSequence.fasta"
#ReadInputFile="22677cDNASingleRead.fastq"
ReadInputFile="22677cDNA.fastq"
AlignmentPrefix="AttemptingGmap"
OutputSubdir=$AlignmentPrefix"/"


#gmap_build -d FullGeneSequence -D /home/minion/MinIONData/cDNA_C04_GMAP $FullGeneReferenceInputFile
#gmap_build -d ExonOnly -D /home/minion/MinIONData/cDNA_C04_GMAP $ExonsOnlyReferenceInputFile

#gmap -f sampe -d FullGeneSequence  -D /home/minion/MinIONData/cDNA_C04_GMAP -4 $ReadInputFile > FullSequenceAlignments.txt
#gmap -f sampe -d ExonOnly  -D /home/minion/MinIONData/cDNA_C04_GMAP -4 $ReadInputFile > ExonsOnlyAlignments.txt

#--idir=/home/minion/MinIONData/TestReadSplit/pass --odir=/home/minion/MinIONData/TestReadSplit/read_extracts -b /home/ben/MUMCScripts/ExtractAllReads/BarcodeSampleMap.txt

