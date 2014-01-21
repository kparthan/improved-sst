import sys, codecs, os, re
if sys.stdout.encoding is None:
        sys.stdout = codecs.open('/dev/stdout', 'w', 'utf-8')

fw = open('experiments-part3.sh','w')

fw.write('STARTM=`date -u "+%s"`\n')
fw.write('line_number=1\n')

MAX_COMPONENTS = 500
cmd = './main-part3 --directory ./spherical_system/class-a/ --mixture --k '
for k in range(2,MAX_COMPONENTS+1):
  #if k % 4 == 0: # part-0
  #if k % 4 == 1: # part-1
  #if k % 4 == 2: # part-2
  if k % 4 == 3: # part-3
    current = cmd + str(k)
    fw.write(current+'\n')
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

os.system('chmod 755 experiments-part3.sh')

