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

import math


# A place to store information about a mapped exon.
class AlignedExon:
    def __init__(self):
        self.referenceSequence=''
        self.alignmentSequence=''
        self.readSequence=''
        self.exonName='UnknownExon'
        self.referenceStartPosition=-1
        self.referenceEndPosition=-1
        
    def assignExonName(self,exonMap):
        # exonMap contains tuples:
        # (ExonName,BeginIndex,EndIndex)
        # I will compare this exon begin/end position with the map file to see if i can name this exon.
        
        specifyPartialExons = False
        
        for exonMapName,exonMapStart,exonMapEnd in exonMap:
           
            # TODO: Shorter exon names.  
            # TODO: Maybe instead of -Partial, I should use '-' and '+' to
            # determine if the mapped sequence is longer or shorter than the reference.
            if closeEnoughEquals(self.referenceStartPosition, exonMapStart):
                if closeEnoughEquals(self.referenceEndPosition, exonMapEnd):
                    #BOTH Start position and end position match.  Ideal.
                    self.exonName = exonMapName
                else:
                    self.exonName = exonMapName
                    if(specifyPartialExons):
                        self.exonName = self.exonName + '-Partial'
                    #ONLY Start positions match.           
                #Stop looking.    
                break 
            elif closeEnoughEquals(self.referenceEndPosition, exonMapEnd):
                #ONLY End positions match.       
                self.exonName = exonMapName
                if(specifyPartialExons):
                    self.exonName = self.exonName + '-Partial' 
                #Stop looking.
                break     
            else:
                #Nothing matches.  Carry on.
                pass
            

            
        

def closeEnoughEquals(value1,value2):
    # This method returns true if value1 and value2 are close to eachother:
    # within "boundaryComparitorBuffer" of eachother    
    #Sometimes the indices are off by a few bases, i need some wiggle room.  5 bases off is okay? 10?
    # TODO: Somethign better with this boundary comparitor buffer.  Maybe it's a CL parameter.
    boundaryComparitorBuffer = 10    
    return (math.fabs(float(value1)-float(value2)) < boundaryComparitorBuffer)      
        
