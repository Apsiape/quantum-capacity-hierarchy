import subprocess,sys,time,json,psutil
from pathlib import Path
OUT=Path(__file__).resolve().parent
args=sys.argv[1:]; tag=args[0]; command=args[1:]; start=time.monotonic(); peak=0; reason=None
with (OUT/(tag+'.log')).open('w') as log:
 p=subprocess.Popen(command,cwd=OUT,stdout=log,stderr=subprocess.STDOUT)
 proc=psutil.Process(p.pid)
 while p.poll() is None:
  try:
   memory=proc.memory_info().rss+sum(c.memory_info().rss for c in proc.children(recursive=True)); peak=max(peak,memory)
  except psutil.Error: pass
  if peak>1400000000: reason='1.4 GB cap'
  if time.monotonic()-start>600: reason='600 second cap'
  if reason:
   try:
    for c in proc.children(recursive=True): c.kill()
    proc.kill()
   except psutil.Error: pass
   break
  time.sleep(.1)
 code=p.wait()
result={'command':command,'exit':code,'seconds':time.monotonic()-start,'peak_tree_RSS':peak,'killed':reason}
(OUT/(tag+'-resources.json')).write_text(json.dumps(result,indent=2)); print(json.dumps(result)); print((OUT/(tag+'.log')).read_text())
sys.exit(code)
