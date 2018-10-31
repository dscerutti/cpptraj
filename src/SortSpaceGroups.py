import sys, os, numpy

finp = open('SpaceGroups.txt', 'r')
fmem = [ line for line in finp ]
ngrp = 0
nops = 0
grpnames = []
symopers = []
symbound = []

for i, ln in enumerate(fmem):

  # Check if the line signals the start of a new space group
  ttli = ln.split()
  if (len(ttli) == 6 and ttli[0] == "Number" and ttli[1] == "of" and ttli[2] == "Symmetry" and
      ttli[3] == "Operators" and ttli[4] == "="):
    j = i
    keepgoing = True
    ngops = 0
    while (j < len(fmem) and keepgoing):
      ttlj = fmem[j].split()
      if (len(ttlj) > 0 and ttlj[0] == "asymm="):
        keepgoing = False
      if (len(ttlj) > 0 and ttlj[0] == "symmetry="):
        thisop = ttlj[1].split(',')
        if (len(thisop) != 3):
          print "Error, I see: ", thisop
          quit
        T = numpy.tile(0.0, 12)
        for k in range(3):

          # Split the string into numerical and axis components, with signs for each
          thisListOp = []          
          for lttr in thisop[k]:
            if (lttr == '+' or lttr == '-'):
              thisListOp.append(' ' + lttr)
            else:
              thisListOp.append(lttr)
          tlv = ''.join(thisListOp)
          newOp = tlv.split()
          for item in newOp:
            if (item == "X" or item == "+X"):
              T[4*k] = 1.0
            elif (item == "-X"):
              T[4*k] = -1.0
            elif (item == "Y" or item == "+Y"):
              T[4*k + 1] = 1.0
            elif (item == "-Y"):
              T[4*k + 1] = -1.0
            elif (item == "Z" or item == "+Z"):
              T[4*k + 2] = 1.0
            elif (item == "-Z"):
              T[4*k + 2] = -1.0
            elif ('X' in item):
              prod = item.split('X')
              if ('/' in prod[0]):
                frac = prod[0].split('/')
                T[4*k] = float(frac[0]) / float(frac[1])
              else:
                T[4*k] = float(prod[0])
            elif ('Y' in item):
              prod = item.split('Y')
              if ('/' in prod[0]):
                frac = prod[0].split('/')
                T[4*k + 1] = float(frac[0]) / float(frac[1])
              else:
                T[4*k + 1] = float(prod[0])
            elif ('Z' in item):
              prod = item.split('Z')
              if ('/' in prod[0]):
                frac = prod[0].split('/')
                T[4*k + 2] = float(frac[0]) / float(frac[1])
              else:
                T[4*k + 2] = float(prod[0])
            elif ('/' in item and 'X' not in item and 'Y' not in item and 'Z' not in item):
              frac = item.split('/')
              T[4*k + 3] = float(frac[0]) / float(frac[1])
            else:
              print "Error, unrecognized operation %s" % item
        symopers.append(T)
        ngops = ngops + 1
      j = j + 1
    if (ngops != int(ttli[5])):
      print "Error, I see %d ops, contrasted with %d I wanted." % (ngops, int(ttli[5]))
      quit
    grpnames.append(fmem[i-1].split()[1])
    symbound.append((nops, nops + ngops))
    ngrp = ngrp + 1
    nops = nops + ngops
    
    # CHECK
#    print "Symmetry operations for %s" % grpnames[ngrp-1]
#    for j in range(symbound[ngrp-1][0], symbound[ngrp-1][1]):
#      for k in range(3):
#        print "  %9.4f %9.4f %9.4f   %9.4f" % (symopers[j][4*k + 0], symopers[j][4*k + 1],
#                                               symopers[j][4*k + 2], symopers[j][4*k + 3])
#      print
    # END CHECK
    
# Write out C++ code to hard-wire all of these space groups and symmetry operations
print "//" + ("-" * 93)
print "// LoadSpaceGroupSymOps: function for creating a list of symmetry operations based",
print "on the space"
print "//                       group specified in an xtalsymm command.  This code was",
print "created by a"
print "//                       python script and should not be edited by hand."
print "//" + ("-" * 93)
print "Action::RetType Action_XtalSymm::LoadSpaceGroupSymOps()"
print "{"
print "  int i, j, k;"
print
print "  // Allocate space to hold all of the symmetry operations."
print "  // No space group has more than 96."
print "  R = new Matrix_3x3[96 * nCopyA * nCopyB * nCopyC];"
print "  T = new Vec3[96 * nCopyA * nCopyB * nCopyC];"
print "  RefT = new Vec3[96 * nCopyA * nCopyB * nCopyC];"
print "  double dnA = (double)nCopyA;"
print "  double dnB = (double)nCopyB;"
print "  double dnC = (double)nCopyC;"
print "  double Rdata[9];"
print
print "  // Loop over all unit cell replicas, then all possible space groups to find a match."
print "  // The giant branch over space groups goes in the nest loops to cut down on the"
print "  // sheer amount of code.  Writing the same loops for every case would make the file"
print "  // even longer, even though the computer will now have to evaluate the switch for"
print "  // each replica of the unit cell during setup."
print "  for (i = 0; i < nCopyA; i++) {"
print "    double di = (double)i;"
print "    for (j = 0; j < nCopyB; j++) {"
print "      double dj = (double)j;"
print "      for (k = 0; k < nCopyC; k++) {"
print "        double dk = (double)k;"
print "        int iuc = (i*nCopyB + j)*nCopyC + k;"
for i, item in enumerate(grpnames):
  if (i == 0):
    print "        if (spaceGrp.compare(\"%s\") == 0) {" % item
  else:
    print "        else if (spaceGrp.compare(\"%s\") == 0) {" % item
  print "          nops = %d * nCopyA * nCopyB * nCopyC;" % (symbound[i][1] - symbound[i][0]);
  print "          iuc *= %d;" % (symbound[i][1] - symbound[i][0])
  for j in range(symbound[i][0], symbound[i][1]):
    print "          Rdata[0] = %14.10f;  Rdata[1] = %14.10f;  Rdata[2] = %14.10f;" % \
      (symopers[j][0], symopers[j][1], symopers[j][2])
    print "          Rdata[3] = %14.10f;  Rdata[4] = %14.10f;  Rdata[5] = %14.10f;" % \
      (symopers[j][4], symopers[j][5], symopers[j][6])
    print "          Rdata[6] = %14.10f;  Rdata[7] = %14.10f;  Rdata[8] = %14.10f;" % \
      (symopers[j][8], symopers[j][9], symopers[j][10])
    print "          R[%3d + iuc] = Matrix_3x3(Rdata);" % (j - symbound[i][0])
    print "          T[%3d + iuc] = Vec3((%14.10f + di) / dnA, (%14.10f + dj) / dnB, " % \
      (j - symbound[i][0], symopers[j][3], symopers[j][7])
    print "                              (%14.10f + dk) / dnC);" % symopers[j][11]
  print "        }"
print "        else {"
print "          return Action::ERR;"
print "        }"
print "      }"
print "    }"
print "  }"
print "  return Action::OK;"
print "}"
    
