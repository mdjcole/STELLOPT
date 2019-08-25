
import sys
import subprocess

filename = sys.argv[1]
print('processing '+filename)

old_list = ['[math](math)', r'\\ '[:-1], r'\_ '[:-1], r'\^', r'\|\| \|\|', r'\|', r"\`", r"\'", r'\$']
new_list = ['$$', r'\ '[:-1], r"_ "[:-1], r"^", r'| \n |', '|', "`", "`", '$$']

head = ['{% include head.html %} \n']

# get table of contents
toc = subprocess.Popen("../gh-md-toc "+filename, shell=True, stdout=subprocess.PIPE).stdout.read()
toc = toc.split(b'\n') # split into list
toc.pop(2)
toc.insert(1, b'---')
toc[-2] = b'---'
toc = [line.decode('utf8').replace('         ', '   ')+' \n' for line in toc]



with open(filename, 'r') as f:
    contents = f.readlines()
    for i, line in enumerate(contents):
        if '====' in line: # insert toc
            if 'Table of Contents' not in contents[i+3]:
                contents[i+1:i+1] = toc
        for old, new in zip(old_list, new_list):           
            contents[i] = contents[i].replace(old, new)

if 'include head.html' not in contents[0]:
    contents = head + contents
            
with open(filename, 'w') as f:
    f.writelines(contents)

exit
