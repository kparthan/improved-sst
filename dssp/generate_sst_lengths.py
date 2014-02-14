import sys, codecs, os, re
if sys.stdout.encoding is None:
        sys.stdout = codecs.open('/dev/stdout', 'w', 'utf-8')

all_dssp_files = os.listdir('./dssp_assignments/')
fw = open('generate_sst_lengths.sh','w')

fw.write('STARTM=`date -u "+%s"`\n')
fw.write('line_number=1\n')

for i in range(len(all_dssp_files)):
  dssp_file = all_dssp_files[i]
  cmd = './sst_lengths --dssp ./dssp/dssp_assignments/' + dssp_file 
  fw.write(cmd+'\n')
  fw.write('echo $line_number\n')
  fw.write('line_number=$((line_number+1))\n')

fw.write('STOPM=`date -u "+%s"`\n')
fw.write('RUNTIMEM=`expr $STOPM - $STARTM`\n')
fw.write('if (($RUNTIMEM>59)); then\n')
fw.write('TTIMEM=`printf "%dm%ds\\n" $((RUNTIMEM/60)) $((RUNTIMEM%60))`\n')
fw.write('else\n')
fw.write('TTIMEM=`printf "%ds\\n" $((RUNTIMEM))`\n')
fw.write('fi\n\n')
fw.write('echo "Executing "script function" took: $TTIMEM"\n')
fw.close()

os.system('chmod 755 generate_sst_lengths.sh')

