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

import sys
import getopt
import os
import subprocess
from subprocess import Popen, PIPE, STDOUT
from os.path import join, split
from Bio import SeqIO
import pysam
import StringIO
from StringIO import StringIO

import numpy
#from numpy import random
from matplotlib import pylab
#import matplotlib.pyplot


import SplicingResult
import AlignedExon


SoftwareVersion = "Alternative-Splicing-Analyzer Version 1.0"

def usage():
    print("usage:\n" + 
    "\tThis description is clearly not correct.  This needs to be fixed.\n\n" + 

    "\tOptions:\n" +  
    "\t-i\t--idir   \tInput Directory (required)\n" +  
   
    "\n\tSee README.MD for instructions on how to set up an anaconda environment for this script\n"
    )    

def loadInputRecords(recordFileName):
    if (".fasta" == recordFileName[-6:] or ".fa" == recordFileName[-3:]):
        FileOutputFormat = "fasta"
    elif (".fastq"== recordFileName[-6:] or ".fq" == recordFileName[-3:]):
        FileOutputFormat = "fastq"
    else:
        FileOutputFormat = "UnknownFormat"
        
    parsedInputReads = SeqIO.parse(recordFileName, FileOutputFormat)    
    alignedSequences = enumerate(parsedInputReads)
    
    return alignedSequences

def readExonMap(exonMapFileName):
    # The exon map is a list of tuples. 
    # A tuple contains
    # (ExonName,StartIndex,EndIndex)
    # The indices are 1-based positions on the full-gene reference sequence,
    # To match IMGT/HLA references. 
    exonMap = []
    
    print ('I will open an exon mapping file:' + exonMapFileName)
    
    try:
        inputFileObject  = open(exonMapFileName, 'r')
        
        for line in inputFileObject.readlines():

            if(len(line.strip()) < 1):
                #print('Empty line:' + line)
                pass            
            # Comments begin with a #
            elif(line.strip()[0] == '#'):
                #print('Commentline:' + line)
                pass
            else:
                #print(line)
                lineSplit = line.strip().split(';')
                exonName = lineSplit[0]
                exonStart = lineSplit[1]
                exonEnd = lineSplit[2]
                exonMap.append((exonName,exonStart,exonEnd))
        
        inputFileObject.close() 
    
    except Exception:
        print 'Problem when reading exon reference.\n'
        print sys.exc_info()[1]
        raise


    return exonMap

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not os.path.isdir(tempDir):
        os.mkdir(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput


# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
    global inputReferenceFilename
    global inputReadFilename
    global outputResultDirectory
    global referenceSequence
    global exonMap
    #global minimumReadLength
    #global maximumReadLength
    
    inputReferenceFilename         = None
    inputReadFilename              = None
    outputResultDirectory          = None
    referenceSequence              = None
    exonMap                        = []
    #minimumReadLength              = 0
    #maximumReadLength              = 0


    # getopt.getopt(..) is a function for parsing args the way smart people do it
    # For More info, use google or 
    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt.getopt(sys.argv[1:]
            ,"hvR:r:o:e:"
            ,["help", "version", "reads=", "reference=" ,"odir=","exonmap="])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                print (SoftwareVersion)
                return False

            #elif opt in ("-i", "--idir"):
            #    inputRootDirectory = arg
            elif opt in ("-o", "--odir"):
                outputResultDirectory = arg
            elif opt in ("-R", "--reads"):
                inputReadFilename = arg
                
            elif opt in ("-e", "--exonmap"):
                exonMap = readExonMap(arg)
                
            elif opt in ("-r", "--reference"):
                inputReferenceFilename = arg
                #Store the reference sequence, we'll need it later.
                for id, record in loadInputRecords(inputReferenceFilename):
                    referenceSequence = str(record.seq)
                #print('Passed Reference Sequence = ' + str(referenceSequence))
            #elif opt in ("-m", "--minlen"):
            #    minimumReadLength = int(arg)
            #elif opt in ("-M", "--maxlen"):
            #   maximumReadLength = int(arg)
            #elif opt in ("-b", "--barcodes"):
            #    barcodeSampleMapFilename = arg
            #    barcodeSampleMap = readBarcodeSampleNames(barcodeSampleMapFilename)

        if(len(sys.argv) < 3):
            print ('I don\'t think you have enough arguments.\n')
            usage()
            return False    

    except getopt.GetoptError, errorMessage:
        print ('Something seems wrong with your commandline parameters.')
        print (errorMessage)
        usage()
        return False

    #print('Input Directory:' + inputRootDirectory)
    print('Output Directory:' + outputResultDirectory)
    #print('Minimum Read Length:' + str(minimumReadLength))
    #print('Maximum Read Length:' + str(maximumReadLength))

    # Quick sanity check.
    #if(len(inputRootDirectory) < 4):
    #    print('Input directory is too short:' + str(inputRootDirectory))
    #    return False
    if(len(outputResultDirectory) < 4):
        print('Output directory is too short:' + str(outputResultDirectory))
        return False
    #if(minimumReadLength < 0):
    #    print('Minimum Read Length should be >= 0 :' + str(minimumReadLength))
    #    return False
    # Actually maximumReadLength can be 0.  That's fine.  No max Length in that case  
    #if(maximumReadLength < 0):
    #    print('Maximum Read Length should be >= 0 :' + str(maximumReadLength))
    #    return False

    return True

    
def buildGmapDatabase():
    print('Building a GMAP Database, from the file:' + inputReferenceFilename)

    try:        
        # Make sure that the  output directory exists.
        if not os.path.exists(outputResultDirectory):
            os.makedirs(outputResultDirectory)
                    
        #Trying to run this command:
        #gmap_build -d FullGeneSequence -D /home/minion/MinIONData/cDNA_C04_GMAP $FullGeneReferenceInputFile
        #how to call a method in python:
        #subprocess.call(["command1", "arg1", "arg2"])
        executionResults = subprocess.call([
            'gmap_build'
            ,'-d','SplicingDatabase'  
            ,'-D',outputResultDirectory
            ,inputReferenceFilename      
            ])
        
        if(executionResults == 0):
            print('Successfully built the gmap splicing database.')
            
        else:
            # I don't know what to do if i can't build a database, lets just crash for now.
            raise Exception('gmap_build was unsuccessful.  \nInput:' + str(inputReferenceFilename) + '\nOutput:' + str(outputResultDirectory))
        


    except Exception:
        print 'Problem when building gmap database.\n'
        print sys.exc_info()[1]
        raise


    
def analyzeAlternativeSplicing():
    # This method will open a fasta or fastq file.
    # Each sequence entry is isolated and compared against the reference sequence.    
    print('Running alternative splicing analysis.  Reads input file:\n' + inputReadFilename)
    splicingResults = []
        
    try:        
        # Make sure that the  output directory exists.
        if not os.path.exists(outputResultDirectory):
            os.makedirs(outputResultDirectory)
   
        # If the read file is a fasta
        if(inputReadFilename.endswith('.fa') or inputReadFilename.endswith('.fasta')):
            readParser = SeqIO.parse(inputReadFilename, "fasta") 
            print ('Reads are in fasta format.')
            # make read parser in fasta format
            
        # If it is fastq
        elif(inputReadFilename.endswith('.fq') or inputReadFilename.endswith('.fastq')):
            readParser = SeqIO.parse(inputReadFilename, "fastq") 
            print ('Reads are in fastq format.')
            # make read parser in fastq format
            
        else:
            readParser = None
            print ('I expected a read file with one of these extensions: .fasta, .fa, .fastq, .fq' )
            raise Exception('Read input invalid format.  Should be fasta or fastq.' )


        for record in readParser:     
            # Pass the read to an analyzeRead method, and store the splicing result.
            splicingResults.append(analyzeRead(record))

        # At this point i have a bunch of splicing results.  I should print it to a output file. 
        spliceSummaryOutputFile = createOutputFile(join(outputResultDirectory,'SplicingSummary.txt'))
        splicedAlignmentOutputFile = createOutputFile(join(outputResultDirectory,'SplicingAlignments.txt'))
        
        
        spliceSummaryOutputFile.write('Splicing Results\n\n')
        
        totalReadCount = len(splicingResults)
        #forwardReadCount = 0
        #unmappedReads = 0
        spliceSummaryOutputFile.write(str(totalReadCount) + ' cDNA Reads Analyzed\n')
        
        #This dictionary stores how many reads have how many exons
        exonCountDictionary={}
        
        #This dictionary stores the exon pattern from each read and counts how many
        exonPatternDictionary={}
        
        # This loop is for counting things.
        for splicingResult in splicingResults:
            #if splicingResult.isForwardStrand:
            #    forwardReadCount += 1
            
            #if(splicingResult.)
                
            exonCount = len(splicingResult.alignedExons)
            if(exonCount in exonCountDictionary):
                exonCountDictionary[exonCount] += 1
            else:
                exonCountDictionary[exonCount] = 1
                
            exonPattern = splicingResult.exonProfile()
            if(exonPattern in exonPatternDictionary):
                exonPatternDictionary[exonPattern] += 1
            else:
                exonPatternDictionary[exonPattern] = 1
            #spliceSummaryOutputFile.write('Hello\n')
            
        #spliceSummaryOutputFile.write(str(forwardReadCount) + ' reads mapped in the forward direction.\n')
        
        
        # Chose 50 because why would there be more exons than that? 
        # TODO This is arbitrary, make this loop better.
        #for i in range(0,50):
        #    if i in exonCountDictionary:
        #        spliceSummaryOutputFile.write(str(exonCountDictionary[i]) + ' reads have ' + str(i) + ' exons.\n')
        
        #Convert the dictionaries to lists and sort em.
        
        exonCountList = [(numExons, exonReadCount) for numExons, exonReadCount in exonCountDictionary.iteritems()]
        exonCountList = sorted(exonCountList, key=lambda tup: tup[1], reverse=True)

        spliceSummaryOutputFile.write(str(exonCountDictionary[0]) + ' reads did not align to the reference.\n\n')
        
        for numExons, exonReadCount in exonCountList:
            if(int(numExons) != 0):
                spliceSummaryOutputFile.write(str(exonReadCount) 
                    + ' reads (' + str( 100.0 * float(exonReadCount) / float(totalReadCount) ) + ' %)'
                    + ' have ' + str(numExons) + ' exons.\n')

        spliceSummaryOutputFile.write('\n')
        
        exonPatternList = [(exonPattern, exonReadCount) for exonPattern, exonReadCount in exonPatternDictionary.iteritems()]
        exonPatternList = sorted(exonPatternList, key=lambda tup: tup[1], reverse=True)
        
        spliceSummaryOutputFile.write('Alternative Splicing Expression Patterns:\n\n')
        
        for exonPattern, exonReadCount in exonPatternList:
            if(len(exonPattern) != 0):
                spliceSummaryOutputFile.write(
                    str(exonPattern) + '\n'
                    + str(exonReadCount) 
                    + ' reads (' + str( 100.0 * float(exonReadCount) / float(totalReadCount) ) + ' %)'
                    + ' follow this pattern.\n\n' 
                    
                    )

        
                
        spliceSummaryOutputFile.write('\n\n')
        
        for splicingResult in splicingResults:
            splicedAlignmentOutputFile.write('>' + splicingResult.id + '\n')
            if splicingResult.isForwardStrand:
                splicedAlignmentOutputFile.write('Forward Read\n')
            else:
                splicedAlignmentOutputFile.write('Reverse Read\n')
                
            exonCount = len(splicingResult.alignedExons)
            if(exonCount > 0):
                splicedAlignmentOutputFile.write('I found ' + str(exonCount) + ' exons in this read:\n')
                splicedAlignmentOutputFile.write(str(splicingResult.exonProfile()) + '\n\n')
                
                #for exonIndex, exon in enumerate(splicingResult.alignedExons):
                                      
                #    splicedAlignmentOutputFile.write(str(exon.exonName) +
                #        ' @ (' + str(exon.referenceStartPosition) + ':' + str(exon.referenceEndPosition) + ')' + 
                #        ': ' + str(len(exon.readSequence)) + ' bp\n')
                   
                #splicedAlignmentOutputFile.write('\n')
                    
                splicedAlignmentOutputFile.write(splicingResult.alignmentText()+ '\n') 
            else:
                splicedAlignmentOutputFile.write('I found 0 exons in this read.\n\n')
              
        spliceSummaryOutputFile.close()
        splicedAlignmentOutputFile.close
        
        #TODO: Make histograms
        histogramListValues = []
        # Make some values for a histogram
        if(len(exonPatternList) > 10):
            patternCount = 10
        else:
            patternCount = len(exonPatternList)
        #    
        for patternIndex in range(0,patternCount):
            #for countIndex in range(0,int(exonPatternList[patternIndex][1])):
            
                histogramListValues.append(exonPatternList[patternIndex])
        #histogramListValues = exonPatternList[0:patternCount]
        
        for histListo in histogramListValues:
            print('list of len ' + str(len(histogramListValues)) + ' has this value:' + str(histListo))
            

        createHistograms(histogramListValues)

    except Exception:
        print 'Problem when running alternative splicing analysis.\n'
        print sys.exc_info()[1]
        raise
 
def analyzeRead(record):    
    # This method takes a single fasta or fastq record, 
    # and aligns the read against the reference sequence using gmap.
    # The gmap output is sent to a results processing method.
    try:
        currentSplicingResult = SplicingResult.SplicingResult()
        currentSplicingResult.id = record.id
        currentSplicingResult.sequence = str(record.seq)
        print('\n\n\nAnalyzing this read:' + str(record.id))
        #print('It has this sequence:' + str(record.seq))
        
        
        AlignmentOutputFolder = join(outputResultDirectory,'Alignments')        
        if not os.path.isdir(AlignmentOutputFolder):
            os.mkdir(AlignmentOutputFolder)

        # Shouldn't need to do this.  Piping from GMAP to pysam is a huge pain.  
        # So I'll write a sam file and then read it, then delete it.
        # Nah, why not just save it?
        samFileName = join(AlignmentOutputFolder,str(record.id) + '.sam')
        temporarySamAlignmentFile = open(samFileName,'w')

        # Open a process that will run gmap.
        p = Popen(
            [
                'gmap'
                ,'-f','sampe' # sampe is sam format, with flags for paired end reads or something.  I don't know.
                ,'-d','SplicingDatabase'  
                ,'-D',outputResultDirectory
                #,'-4' # -4 means show each exon alignment on a different line.  This is human readable.  But I want a sam file.
            ]
            , stdout=temporarySamAlignmentFile, stdin=PIPE, stderr=PIPE
        )
        
        # Prepare the record as a fasta, and pipe the text into the gmap process.
        currentFastaString = record.format('fasta')

        # Input=fasta read.  Output=temporary sam file      
        stdout_data, stderr_data = p.communicate(input=currentFastaString)
        temporarySamAlignmentFile.close()
        
        print('Done with GMAP alignment, next I will open the alignment file.')

        # Create the pysam object which we can work with. 
        # Give this a name, and put the alignment in the output directory. 
        #currentAlignment = pysam.AlignmentFile(samFileName, 'r')
        
        # Analyze splicing.
        currentSplicingResult.processGmapInput(samFileName, referenceSequence, exonMap)
        
        #currentAlignment.close()
        # I'm nervous about deleting files for some reason.  Hope this doesn't mess something up.
        # Don't do this.  We should just keep the alignment files.  They might be interesting, and don't take up much space.
        #os.remove('tempsamalignment.sam')

        return currentSplicingResult
    
    
    except Exception:
        print ('Problem when processing a read:' + str(record.id) + '\n')
        print sys.exc_info()[1]
        raise
  
  
  
  
  
def createHistograms(values):
    print('Let us now create some graphs')
    print('Values:\n')
    print(str(values))
    
    #numberCategories = len(values)
    
    for value in values:
        expressionPattern, patternCount = value
        if(len(expressionPattern.strip()) < 1):
            #values[barIndex] = ('No Mapped Exons', patternCount)
            # TODO: instead of renaming it No Mapped Exons, maybe I should just remove this value.
            print('THE VALUE OF VALUES IS:' + str(values))
            print('THE TYPE of VALUES IS:' + str(type(values)))
            # IT isa list, i should be able to just remove the value.  Do it.
            values.remove(value)
        
    
    #patternValues = zip(*values)[0]
    #patternCounts = [float(i) for i in zip(*values)[1]]
    
    countBars = len(values)

    #means_men = (20, 35, 30, 35, 27)
    #means_men = patternCounts
    #std_men = (2, 3, 4, 1, 2)
    
    #means_women = (25, 32, 34, 20, 25)
    #std_women = (3, 5, 2, 3, 3)
    
    colors = ['b','m','g','c','y','r','k','gold','springgreen','saddlebrown']
    
    fig, ax = pylab.subplots()
    
    indices = numpy.arange(countBars)
    bar_width = .9
    
    opacity = 0.9
    error_config = {'ecolor': '0.3'}
    
    
    
    #def labelBar(rect):
        #for rect in rects:
    #    height = rect.get_height()
        #ax.text('LABEL', ha='center', va='bottom')
    
    #for exonBoundary in exonBoundaries:
    #    exonHandle, = pylab.plot(exonBoundary, exonheight, 'r', label='Exons', linewidth=2.0)
    barHandles = []
    
    
     
    for barIndex in indices:
        expressionPattern, patternCount = values[barIndex]
        #means_men = patternCounts
        
        #if(len(expressionPattern.strip()) < 1):
            #values[barIndex] = ('No Mapped Exons', patternCount)
            # TODO: instead of renaming it No Mapped Exons, maybe I should just remove this value.
         #   print('THE VALUE OF VALUES IS:' + str(values))
         #   print('THE TYPE of VALUES IS:' + str(type(values)))
            # IT isa list, i should be able to just remove the value.  Do it.
            #a.remove(values[barIndex])
            
        
        #print ('this bar index is ' + str(barIndex))
        #print ('for this bar, the pattern is ' + str(expressionPattern))
        #print ('the count is ' + str(patternCount))
    
        barHandle, = pylab.bar([barIndex], [patternCount], bar_width,
            alpha=opacity,
            #color='b',
            color=colors[barIndex],
            #yerr=std_men,
            error_kw=error_config,
            label=expressionPattern)
        
        #labelBar(barHandle)
        
        # Give the bar a label.
        #TODO: I want to give a percentage label to each bar.  It's tough to do.
        # https://matplotlib.org/examples/api/barchart_demo.html
        

        
        barHandles.append(barHandle)
    


    pylab.xlabel('Expression Patterns')
    pylab.ylabel('Read Count')
    pylab.title('Expression Analysis')
    #pylab.xticks(indices + bar_width / 2, ('A', 'B', 'C', 'D', 'E','F','G','H','I','J'))
    #pylab.xticks = indices + 1
    
    patternValues = zip(*values)[0]
    pylab.legend(handles=barHandles, labels=patternValues )
    #pylab.legend()
    
    pylab.tight_layout()
    pylab.show()


if __name__=='__main__':
    try:        
        if(readArgs()):
            print('Commandline arguments look fine.  Starting to analyze splice sites.')
            
            #TODO: Maybe a commandline parameter to skip creating the database?
            #buildGmapDatabase()
            analyzeAlternativeSplicing()
            print('\nDone.  Ben is the greatest programmer.\n')
    
        else:
            print('\nI give up, something\'s not right with the args.') 


    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything except crying
        print 'Unexpected problem during execution:'
        print sys.exc_info()[1]
        raise


