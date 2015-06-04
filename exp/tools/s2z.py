import os
import sys

def freplace(fname, from_pattern, to_pattern):
  data = open(fname).read().replace(from_pattern, to_pattern)
  open(fname, "w").write(data)

def main(argv):
  for fname in argv[1:]:
    head, rest = os.path.split(fname)
    zfname = os.path.join(head, "z" + rest[1:]) # assumes "z" should go in the beginning of file name

    freplace(fname, "@date September 2014
    os.system("python %s -p z --file %s" % (os.environ.get("PATH_TO_CODEGEN_PY", "codegen.py"), fname))
    freplace(zfname, "@generated z", "@precisions normal z -> s d c")

  return 0

sys.exit(main(sys.argv))
