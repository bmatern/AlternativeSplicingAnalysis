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

import pysam
import sys
from AlignedExon import AlignedExon

# Information about a MinION Read's splice variant.
class SplicingResult:
    def __init__(self):
        self.id=''
        self.alignedReferenceSequence=''
        self.alignedReadSequence=''
        self.isForwardStrand=True
        self.alignedExons=[]

        
    def processGmapInput(self, alignmentSamFilename, referenceSequence, exonMap):
        try:
            #self.sequence = 'AGCT'
            print('Processing samfile input.  Read=' + str(self.id))


            print('Step 1 is to convert to a bam file.')
            alignmentBamFilename = alignmentSamFilename.replace('.sam','.bam')
            infile = pysam.AlignmentFile(alignmentSamFilename, 'r')
            outfile = pysam.AlignmentFile(alignmentBamFilename, "wb", template=infile)
            for s in infile:
                outfile.write(s)
            infile.close()
            outfile.close()
            
            print('Step 2 is to index the bam file.')
            pysam.index(alignmentBamFilename)
            
            
            
            print('Step 3 is to open the bam file and look for sequence alignments.')            
            bamfile = pysam.AlignmentFile(alignmentBamFilename, 'rb')
            
            
            #alignedSequenceIterator = bamfile.fetch('HLA-C*04_01_01_01', 1, 1200)
            alignedSequenceIterator = bamfile.fetch()
            
            
            alignedSequenceCount = 0
            
            # In general there is only a single aligned sequence here.  
            # But the alignment has lots of different positions
            # I must iterate through the data and pull out what is relevant.
            # If there is nothing here, then the read didn't align well.
            for alignedSequence in alignedSequenceIterator:
                alignedSequenceCount += 1
                
                #TODO: this isn't quite right.  Sometimes the program is just not sure. if it is forward and backward. 
                # Perhaps a case statement for forward,reverse,unknown,no alignment
                self.isForwardStrand = not alignedSequence.is_reverse
                
                #print('Inside the iterator')
                print ('Read # ' + str(alignedSequenceCount))
                #print ('Is this read mapped to the forward (sense) strand?:' + str(self.isForwardStrand))


                self.alignedReferenceSequence = referenceSequence[alignedSequence.reference_start:alignedSequence.reference_end]
                #print ('reference alignment sequence:' + self.alignedReferenceSequence)
                
                self.alignedReadSequence=str(alignedSequence.query_alignment_sequence)                
                
                if(self.isForwardStrand):
                    print('SENSE STRAND')
                else:
                    print('ANTISENSE STRAND.')
                       
                print('I will try to parse the cigar tuples.')
                    
                    
                # AG TAGCT
                # ||-| |-|
                # AGCTCG T
                    
                self.alignedExons=[]
                
                
                currentAlignedExon = AlignedExon()

                referenceIndex=0
                readIndex=0
                

                cigarTuples = alignedSequence.cigartuples
                for cigarTuple in cigarTuples:
                    #print ('Tuple:' + str(cigarTuple))
                    
                    if(str(cigarTuple[0]) == '0'):
                        regionLength = int(cigarTuple[1])
                        #print('This tuple is for a Match segment of length ' + str(cigarTuple[1]))
                        
                        # Assign the reference sequence begin for this aligned exon.
                        # alignedSequence.reference_start and referenceIndex are 0-based, but I want the stored reference sequences 1-based.
                        # 1-based means I can compare with IMGT/HLA reference values.  
                        
                        if(currentAlignedExon.referenceStartPosition == -1): # If it is unassigned
                            currentAlignedExon.referenceStartPosition = alignedSequence.reference_start + referenceIndex + 1
                        
                        
                        for i in range(0,regionLength):
                            currentAlignedExon.referenceSequence += self.alignedReferenceSequence[referenceIndex]
                            currentAlignedExon.readSequence += self.alignedReadSequence[readIndex]
                            
                            if(self.alignedReferenceSequence[referenceIndex] == self.alignedReadSequence[readIndex] ):
                                currentAlignedExon.alignmentSequence += '|'
                            else:
                                currentAlignedExon.alignmentSequence += ' '
                            
                            referenceIndex += 1
                            readIndex += 1
                         
                        # At the end of a match section.   
                        # assign the end reference position. 
                        # for converting from 0 to 1 based value, nothing needs to be done for the end index.
                        # This will be overwritten if there are more match sections. That's okay.
                        currentAlignedExon.referenceEndPosition = alignedSequence.reference_start + referenceIndex
                            
          
                        
                    elif(str(cigarTuple[0]) == '1'):
                        #print('This tuple is for a Insertion segment of length ' + str(cigarTuple[1]))
                        
                        regionLength = int(cigarTuple[1])
                        
                        # In a insertion, insert a gap from the reference, and display the character in the read.
                        for i in range(0,regionLength):
                            currentAlignedExon.referenceSequence += ' '
                            currentAlignedExon.readSequence += self.alignedReadSequence[readIndex]
                            currentAlignedExon.alignmentSequence += '-'
                            
                            #referenceIndex += 1
                            readIndex += 1                            
                        
                    elif(str(cigarTuple[0]) == '2'):
                        #print('This tuple is for a Deletion segment of length ' + str(cigarTuple[1]))
                        
                        regionLength = int(cigarTuple[1])
                        
                        # In a deletion, display the character from the reference, and insert a gap in the read.
                        for i in range(0,regionLength):
                            currentAlignedExon.referenceSequence += self.alignedReferenceSequence[referenceIndex]
                            currentAlignedExon.readSequence += ' '
                            currentAlignedExon.alignmentSequence += '-'
                            
                            referenceIndex += 1
                            #readIndex += 1
                        
                    elif(str(cigarTuple[0]) == '3'):
                        #print('This tuple is for a BAM_CREF_SKIP segment of length ' + str(cigarTuple[1]))
                        #print('This is an intron, so we can start a new EXON.')
                        
                        regionLength = int(cigarTuple[1])
                        

                        
                        currentAlignedExon.assignExonName(exonMap)
                        self.alignedExons.append(currentAlignedExon)
                        

                        currentAlignedExon = AlignedExon()
                        
                        #Feed out a length of the reference sequence.  
                        referenceIndex += regionLength
                        
                    elif(str(cigarTuple[0]) == '4'):
                        #print('This tuple is for a BAM_CSOFT_CLIP segment of length ' + str(cigarTuple[1]))
                        #print('Soft clips are unaligned sections before the exon alignments.  Nothing to do.')
                        pass
                    
                    # I don't expect to see these bam flags.  If I do, throw an exception.   Figure out what it means
                    elif(str(cigarTuple[0]) == '5'):
                        print('This tuple is for a BAM_CHARD_CLIP segment of length ' + str(cigarTuple[1]))
                        raise(Exception('The developer forgot to handle this cigar tuple:' + str(cigarTuple[0])))
                    elif(str(cigarTuple[0]) == '6'):
                        print('This tuple is for a BAM_CPAD segment of length ' + str(cigarTuple[1]))
                        raise(Exception('The developer forgot to handle this cigar tuple:' + str(cigarTuple[0])))
                    elif(str(cigarTuple[0]) == '7'):
                        print('This tuple is for a BAM_CEQUAL segment of length ' + str(cigarTuple[1]))
                        raise(Exception('The developer forgot to handle this cigar tuple:' + str(cigarTuple[0])))
                    elif(str(cigarTuple[0]) == '8'):
                        print('This tuple is for a BAM_CDIFF segment of length ' + str(cigarTuple[1]))
                        raise(Exception('The developer forgot to handle this cigar tuple:' + str(cigarTuple[0])))
                    else:
                        print('I have no idea what is going on with this tuple, so I will throw a fit.')
                        raise(Exception('Unknown Tuple identifier:' + str(cigarTuple[0])))
                
                currentAlignedExon.assignExonName(exonMap)
                self.alignedExons.append(currentAlignedExon)
                 
                #print('Here are the alignments:\n')
                #print(self.alignmentText())
                        
                #else:
                   # print('This read is on antisense strand, I will figure this out later.')
             
             
                
            bamfile.close()
            
            if(alignedSequenceCount == 1):
                print ('There was one aligned sequence in that bam file.  This is what I expected.')
            elif(alignedSequenceCount == 0):
                print ('There was 0 aligned sequence in that bam file.  I guess nothing aligned to the reference.')
            else:
                print ('Multiple sequences bound to the same reference in this bam file.  I did not expect this.')
                raise(Exception('Multiple sequences bound to the same reference file.'))
                
        except Exception:
            print('Problem when processing the SAM file:' + self.id)
            print sys.exc_info()#[1]
            
            #Extremely dangerous, we'll lose some information this way,but a bam crashes sometimes and i don't know why.
            #TODO: Something better with this exception.  Right now I'm ignoring it.
            # I should pass-through the corrupt-bam error unless i can find a file.
            #raise
    
    def alignmentText(self):
        
        
                         #   splicedAlignmentOutputFile.write(str(exon.exonName) +
                      #  ' @ (' + str(exon.referenceStartPosition) + ':' + str(exon.referenceEndPosition) + ')' + 
                      #  ': ' + str(len(exon.readSequence)) + ' bp\n')
        
        #Generate some verbose text to show what exons are in this splicing result
        
        
        alignmentText = ''
           
        for exon in self.alignedExons:                      
            alignmentText += (str(exon.exonName) +
                ' @ (' + str(exon.referenceStartPosition) + 
                ':' + str(exon.referenceEndPosition) + ')' + 
                ': ' + str(len(exon.readSequence)) + ' bp\n')
            alignmentText += (exon.referenceSequence + '\n')
            alignmentText += (exon.alignmentSequence + '\n')
            alignmentText += (exon.readSequence + '\n')
            alignmentText += ('\n')
            
        return alignmentText
    
    def exonProfile(self):
        # A short string with a list of exon names.  To be used for categorization, quantification and sorting        
        exonProfile = ''
           
        for exon in self.alignedExons:  
            exonProfile +=  (str(exon.exonName) + ' : ')
            
        #trim the last comma off.
        exonProfile = exonProfile[0:len(exonProfile)-3].strip()
            
        return exonProfile
        
        
        
