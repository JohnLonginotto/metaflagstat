#!/usr/bin/env python
# coding=utf8
# -*- coding: utf8 -*-

###############################
## shallowTools.metaflagstat.py
###############################
## This tool takes n bam/sam files as arguments, and for each file creates a childprocess to read it and count all unique flags.
## These flags are then bitwise compared to the 4094 possible combinations of flags (0 to 111111111111), and the value for every possible combination created (even if the occurence of that flag doesn't exist)
## So say 000000000001 (read is paired) occurs for every single read in a bam, but no read has ONLY 000000000001 set, here the final array will actually contain a value for 000000000001.
## By breakig down/adding up the reads like this, you can very easily add/subtract/divide/multiple flag counts for all possible queries without having to re-read the bam file.
## This data, for every sample, is then added to an HTML file with Highcharts.js, where you can select the flags to include/exclude with AND/OR comparitors. It's how flagstat should have been.
## If you give it a bam file, and have samtools AND pysam installed, it'll use pysam (about 2x faster than samtools view -> python alone).
## If you give it a sam file, it'll just read the file.

samtools = '/path/to/samtools'


## JOHN:
## counts = {} or collections.counter? speedtest.

'''
Flag        Int     Description
0x0001      1       the read is paired in sequencing
0x0002      2       the read is mapped in a proper pair
0x0004      3       the query sequence itself is unmapped
0x0008      4       the mate is unmapped
0x0010      5       strand of the query (1 for reverse)
0x0020      6       strand of the mate
0x0040      7       the read is the first read in a pair
0x0080      8       the read is the second read in a pair
0x0100      9       the alignment is not primary
0x0200      10      the read fails platform/vendor quality checks
0x0400      11      the read is either a PCR or an optical duplicate
0x0800      12      the read is a supplimentary alignment
'''
 
import subprocess
import collections
import sys
import os
import fileinput
import json

def bamCheck(filename):
	xamfile = open(filename, 'rb')
	try:
		grab = 2048
		while 1:
			section = xamfile.read(grab)
			if '\0' in section: return True
			if len(section) < grab: break
	finally: xamfile.close()
	return False

files = sys.argv[1:]
list = {}
x = 0
for file in files:
	fileName = os.path.basename(file)
	x += 1
	if file == '--samtools':
		sys.argv = ''	## dumb, but has to be done so fileinput doesn't take args instead of stdin.
		counter = {}
		for line in fileinput.input():
			try:
				counter[line.split('\t')[1]] += 1
			except KeyError:
				counter[line.split('\t')[1]] = 1
		answers = {}
		for qflag in range(0,4094):
			temp = 0
			for flag, count in counter.items():
				flag = int(flag)
				if (flag & qflag) == qflag:
					temp += int(count)
			if temp != 0: answers[qflag] = temp
		print json.dumps(answers)
		exit()

	if file == '--pysam':
		import pysam
		file = sys.argv[2]
		counter = collections.Counter()
		counter = {}
		sortedBamfile = pysam.Samfile(file, "rb")
		for alignedread in sortedBamfile:
			try:
				counter[alignedread.flag] += 1
			except KeyError:
				counter[alignedread.flag] = 1 ## If we use collections.Counter, we dont need this try/except, nor this line. We can just jam things in and if it's not there it adds it. 
		answers = {}
		for qflag in range(0,4094):
			temp = 0
			for flag, count in counter.items():
				flag = int(flag)
				if (flag & qflag) == qflag:
					temp += int(count)
			if temp != 0: answers[qflag] = temp
		print json.dumps(answers)
		exit()

	## OK, to be here user must have evoked script with just a list of sam/bed files. Find if binary, and apply appropriate method:
	thisFile = os.path.abspath(__file__)
	if bamCheck(file):
		#print 'File is BAM'
		## Assuming file is BAM, trying to import pysam:
		try:
			import pysam
			pysam = True
		except:
			pysam = False
		if pysam:
			#print 'Using pysam..'
			list[fileName] = subprocess.Popen(thisFile + ' --pysam ' + file, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')
		else:
			#print 'Using samtools..'
			list[fileName] = subprocess.Popen(samtools + ' view ' + file + ' | python ' + thisFile, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')
	else:
		#print 'File is SAM'
		list[fileName] = subprocess.Popen('cat ' + file + ' | python ' + thisFile + ' --samtools', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')
		## THIS IS ONLY GOING TO WORK IF SAM FILES ARE THE SAME AS THE OUTPUT OF samtools view ./somebam - check it!

results = []
for file, processOut in list.items():
	result = json.loads(processOut.communicate()[0][:-1])
	results.append((file,result))
	allData = json.dumps(results)

print '''
<html>
	<head>
	<title>metaflagstat</title>
	<meta name="author" content="John Longinotto - longinotto@immunbio.mpg.de">
	<!-- If you want to load the Javascript from the internet and remove it from this file, then just uncomment these lines: -->
	<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.8.2/jquery.min.js"></script>
	<script src="http://code.highcharts.com/highcharts.js"></script> 
	<script src="http://code.highcharts.com/modules/exporting.js"></script>
	<!-- -->

    </script>

	<script type="text/javascript">
	var allData ='''
print allData + '''

	allKeys = {A:1, C:2, E:4, G:8, I:16, K:32, M:64, O:128, Q:256, S:512, U:1024, W:2048, Z:0,
					B:'', D:'', F:'', H:'', J:'', L:'', N:'', P:'', R:'', T:'', V:'', X:''}
	$(function () { 
		$('#container').highcharts({
		credits: {
			text: 'metaflagstat.py - ac.gt <br> Charts made by Highcharts.js (which, I am required to inform you, is not free for commercial use)',
			href: 'http://ac.gt/metaflagstat'
		},
		chart: {
			type: 'bar'
		},
		title: {
			text: 'Read Count Statistics'
		},
		xAxis: {
			categories:'''
sampleList = []
for result in results:
	sampleList.append(result[0])
print sampleList
print '''
		},
		yAxis: {
			title: {
				text: 'Reads'
			}
		},
		series: [ '''
## Create a series array for every QUESTION (bars in bar graph). Data is ordered by sampleList/files.
questions = [('All Reads',0)]
for question in questions:
	questionTitle = question[0]
	questionInt = str(question[1])
	print '{ name: "' + str(questionTitle) + '", data:',
	data = []
	for result in results:
		data.append(int(result[1][questionInt]))
	print data,
	if question == questions[-1]:
		print '''}
'''
	else:
		print '''},
'''

print ''' ]
	});
});

$('#addData').live('click', function () {
    // Parse the question.
    questionName = $('#dataName').val()                 // Value for series name
    protoQuestion = $('#question').val().toUpperCase()  // Get user's query, converts everything to uppercase.
    questions = [];                                     // algebra broken down into question subsets, like "AB+C" would be "intersect of AB", plus "the whole of C".
    x = 0
    questions[x] = undefined;
    for (var i = 0, len = protoQuestion.length; i < len; i++) {         // iterate through every character in user's query
        character = protoQuestion.charAt(i);
        // Some checks..
        if (character == ' ') { continue; }
        else if ('ABCDEFGHIJKLMNOPQRSTUVWXZ0123456789+-/*'.indexOf(character) === -1) { alert('The symbol ' + character + 'confuses me! Please leave it out :('); return }

        else if ('ABCDEFGHIJKLMNOPQRSTUVWXZ'.indexOf(character) !== -1) {
            // Got Algebra!
            if (typeof questions[x] == 'undefined')   { questions[x] = Number(allKeys[character]); }                                    // Previously had nothing.
            else if (typeof questions[x] == 'number') { questions[x] += Number(allKeys[character]); }                                   // Previously had other algebra, so add bitflags.
            else if ('-+/*'.indexOf(questions[x]) !== -1) { x+=1; questions[x] = Number(allKeys[character]); }                          // Previously had operator. Move on and add.
            else if (typeof questions[x] == 'string') { x+=1; questions[x] = '*'; x+=1; questions[x] = Number(allKeys[character]);  }   // Previously had str(int). Move, add *, move and add.
            else { alert('Massive Internal Error. Possibly Fatal. I have no idea what happened. Please tell John what you typed :S'); return }
        } 
        else if ('-+/*'.indexOf(character) !== -1) {
            // Got an operator!
            if (typeof questions[x] == 'undefined')   {                                                                         // Previously had nothing.
                if (character == '-') { questions[x] = '-';}                                                                    // like -A or -1000. This is OK.
                else if (character == '/' || character == '*') { alert('You can not start with * or / !!'); return }   			// like /A or *A. Not cool bro.
            }                                                                                                                   // No need to do anything for +. +A or +100 is a given.
            else if (typeof questions[x] == 'number') { x+=1; questions[x] = character; }                                       // Previously had algebra. move, add.
            else if ('-+/*'.indexOf(questions[x]) !== -1) {                                                                     // Previously had another operator. Do:
                if(questions[x] == '-' && character == '+') { continue }                                                        // -+ = - [not needed]
                else if(questions[x] == '-' && character == '-') { questions[x] = '+' }                                         // -- = +
                else if(questions[x] == '-' && character == '*') { alert('You can not do -*.. '); return }                      //
                else if(questions[x] == '-' && character == '/') { alert('You can not do -/.. '); return }                      //

                else if(questions[x] == '+' && character == '+') { continue }                                                   // ++ = + [not needed]
                else if(questions[x] == '+' && character == '-') { questions[x] = '-' }                                         // +- = -
                else if(questions[x] == '+' && character == '*') { alert('You can not do +*.. '); return }                      //
                else if(questions[x] == '+' && character == '/') { alert('You can not do +/.. '); return }                      //

                else if(questions[x] == '*' && character == '+') { continue }                                                   // [not needed]
                else if(questions[x] == '*' && character == '-') { x+=1; questions[x] = '-'; }                                  // if operator * or /, move one, if -, apply to NEXT then * or /.
                else if(questions[x] == '*' && character == '*') { /*questions[x] = '**'*/ alert('i dont support that yet.'); return}// ** = ** (to the power of) - would have to add ** to indexOf, and add a new block of 4 options for ** (or renamed to ^?, same thing.)
                else if(questions[x] == '*' && character == '/') { alert('You can not do */.. '); return }                      //

                else if(questions[x] == '/' && character == '+') { continue }                                                   // 
                else if(questions[x] == '/' && character == '-') { x+=1; questions[x] = '-'; }                                  // 
                else if(questions[x] == '/' && character == '*') { alert('You can not do /*.. '); return }                      // 
                else if(questions[x] == '/' && character == '/') { alert('You can not do //.. '); return }                      // 
            }
            else if (typeof questions[x] == 'string') { x+=1; questions[x] = character; }                                       // Previously had str(int). move, add. consciously verbose.
            else { alert('Massive Internal Error. Possibly Fatal. I have no idea what happened. Please tell John what you typed :S'); return }
 
        } else if ('0123456789'.indexOf(character) !== -1) {
            // Got a number (expressed as a string)!
            if (typeof questions[x] == 'undefined')   { questions[x] = character; }                                             // Previously had nothing.
            else if (typeof questions[x] == 'number') { x+=1; questions[x] = '*'; x+=1; questions[x] = String(character); }     // Previously algebra, so x+=1, add *, x+=1, add this.
            else if (questions[x].match(/^[0-9]+$/) != null) { questions[x] = String(questions[x]) + String(character); }       // Previously another str(int), so just concat them :)
            else if ('-+/*'.indexOf(questions[x]) !== -1) { x+=1 ; questions[x] = String(character); }                          // Previously had operator. x+=1, add this.
        }
     }
// operator by default 'undefined'. if see a +, just add next to counts. if you see a -, add negative of next to counts. if you see a * or /, THEN change the operator to * or /.
// if you get to a string or int and the operator IS NOT undefined, then it's time to apply operator to counts. then set operator back to undefined.
    questionResults = [];
    j = 0;                                                																		// j += 1 says we're done with this FILE, move on.
    // For every FILE:
    $.each(allData, function( key, file ) {
        operator = undefined;                         																 			// by default, no operator is set, and we add everything.
        fileData = file[1];
        // For every SUBQUESTION:
        for ( var p = 0, len = questions.length; p < len; p++) {
            //alert('subquestion ' + p + ' = ' + questions[p]);
            if (questionResults[j] == undefined) {
                // OK, we're doing are first subquery for this file.
                if(questions[p] == '-') {                                                           // First value is (the only permissible first operator) "-". Jump forwards, 0 - next value.
                    p+=1;
                    if (typeof questions[p] == 'number')   { questionResults[j] = 0-Number(fileData[questions[p]]); }
                    else if (typeof questions[p] == 'string') { questionResults[j] = 0-Number(questions[p]); }
                    else { alert('explosions, etc.'); return }
                }
                else if(typeof questions[p] == 'string') { questionResults[j] = Number(questions[p]); }     					// First value is a str(int).
                else if (typeof questions[p] == 'number') { questionResults[j] = Number(fileData[questions[p]]) }  				// First value is a bitflag. Grab counts and pop in.
                else { alert('more explosions!'); return }
            } else {
                // OK, so there's deffo something in q[j]. Got to play it safe now...
                if (operator != undefined) {                                                                     				// if operator set from last round, apply it!
                    if (questions[p] == '-') {                                                                      			// unless this is a '-', in which case negate then apply...
                        if      (typeof questions[p] == 'string') { tempResult = 0-Number(questions[p]); }
                        else if (typeof questions[p] == 'number') { tempResult = 0-Number(fileData[questions[p]]); }
                        if      (operator == '*') { questionResults[j] = questionResults[j] * tempResult; }
                        else if (operator == '/') { questionResults[j] = questionResults[j] / tempResult; }
                        operator = undefined;                                                                       			// and reset operator.
                    }
                    else {
                        if      (operator == '*') { 
                            if      (typeof questions[p] == 'string') { questionResults[j] = questionResults[j] * Number(questions[p]) }
                            else if (typeof questions[p] == 'number') { questionResults[j] = questionResults[j] * Number(fileData[questions[p]]) }
                        }
                        else if (operator == '/') { 
                            if      (typeof questions[p] == 'string') { questionResults[j] = questionResults[j] / Number(questions[p]) }
                            else if (typeof questions[p] == 'number') { questionResults[j] = questionResults[j] / Number(fileData[questions[p]]) }
                        }
                        operator = undefined;                                                                       			// and reset operator.
                    }
                }
                else if (questions[p] == '*' || questions[p] == '/') { operator = questions[p]; }                   			// Operator not previously set, but is now operator! set it!
                else if (questions[p] == '-') { 
                    p+=1; 
                    if      (typeof questions[p] == 'string') { questionResults[j] -= Number(questions[p]) }
                    else if (typeof questions[p] == 'number') { questionResults[j] -= Number(fileData[questions[p]]) }
                    else    { alert('so many explosions...'); }
                }
                else if (questions[p] == '+') { 
                    p+=1; 
                    if      (typeof questions[p] == 'string') { questionResults[j] += Number(questions[p]) }
                    else if (typeof questions[p] == 'number') { questionResults[j] += Number(fileData[questions[p]]) }
                }
                else if (typeof questions[p] == 'string') { questionResults[j] += Number(questions[p]); }                   	// Got a str(int)! Add to count.
                else if (typeof questions[p] == 'number') { questionResults[j] += Number(fileData[questions[p]]); }
                else {
                    alert('huh?');
                    // Got an int, which means it was originally letters like 'A' or 'DGQ'
                    if (operator == ' ') { questionResults[j] = Number(fileData[questions[p]]); }
                    if (operator == '+') { questionResults[j] += Number(fileData[questions[p]]); }
                    if (operator == '-') { questionResults[j] -= Number(fileData[questions[p]]); }
                    if (operator == '/') { questionResults[j] = questionResults[j] / Number(fileData[questions[p]]); }
                    if (operator == '*') { questionResults[j] = questionResults[j] * Number(fileData[questions[p]]); }
                }
            }
            //alert(questionResults[j]);
        }
        j += 1;
    });
    chart = $('#container').highcharts();
    chart.addSeries({
        name: questionName,
        data: questionResults
    });
});



</script>


</head>
<body>
<div id="questions" style="width:100%; height:10%; border-bottom: 3px solid grey;">
    <input type="submit" id="addData" value="Add Data!">
    <input type="text" id="dataName" value="Data Name">
    <input type="text" id="question" value="question">
</div>
<div id="container" style="position:relative; float:left; width:85%; height:90%;"></div>
<div id="dictionary" style="width:15%;height:90%;text-align:right;float:right;padding-top:2em">
<pre>
KEY:
Read Paired = A
Read Not Paired = B
Read Mapped (Proper Pair) = C
Read Not Mapped (in P.P.) = D
Read Unmapped = E
Read Not Unmapped = F
Mate Unmapped = G
Mate Not Unmapped = H
Read Reverse Strand = I
Read Not Reverse Strand = J
Mate Reverse Strand = K
Mate Not Reverse Strand = L
First In Pair = M
Not First In Pair = N
Second In Pair = O
Not Second In Pair = P
Not Primary Alignment = Q
Primary Alignment = R
Read Fails QC Checks = S
Read Does Not Fail QC = T
Read Is Duplicate = U
Read Is Not Duplicate = V
Supplementary Alignment = W
Not Supplementary = X
ALL READS = Z
</pre>
<div>



  

</body></html>
'''
