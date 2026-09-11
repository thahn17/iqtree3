import re
import sys

scratch = sys.argv[1]

def strip_cache(s):
    s = re.sub(r'(<c:numCache>).*?(</c:numCache>)', r'\1...\2', s, flags=re.DOTALL)
    s = re.sub(r'(<c:strCache>).*?(</c:strCache>)', r'\1...\2', s, flags=re.DOTALL)
    return s

for i in [2, 3, 4]:
    path = f"{scratch}/chartinspect/xl/charts/chart{i}.xml"
    s = open(path, encoding='utf-8').read()
    s = strip_cache(s)
    out = f"{scratch}/chart{i}_skeleton.xml"
    open(out, 'w', encoding='utf-8').write(s)
    print(i, len(s))
