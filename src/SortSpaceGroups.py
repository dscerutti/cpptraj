import sys, os, numpy

finp = open('SpaceGroups.txt', 'r')
fmem = [ line for line in finp ]
ngrp = 0
nops = 0
grpnames = []
symopers = []
symbound = []
grpcomps = []

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

        # Determine the rules for the primary asymmetric unit's geometric bounds
        nrule = 1
        rules = []
        rlstring = ""
        for item in ttlj:
          if (item == "asymm="):
            pass
          elif (item == "and"):
            nrule = nrule + 1
            rules.append(rlstring)
            rlstring = ""
          else:
            rlstring += item
        rules.append(rlstring)
        grpcomps.append([])
        for item in rules:
          tti = item.split("<=")
          if (len(tti) == 2):
            grpcomps[ngrp].append((tti[0], tti[1]));
          elif (len(tti) == 3):
            grpcomps[ngrp].append((tti[0], tti[1]));
            grpcomps[ngrp].append((tti[1], tti[2]));
          else:
            print "Error, there are %d combined comparisons." % len(tti)
            exit(1)

        # CHECK
        print fmem[j].split("\n")[0]
        print grpcomps[ngrp]
        print
        # END CHECK
        
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
  print "          sgID = %d;" % i
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

# Print the preamble to the asymmetric unit definitions function
print "//" + ("-" * 93)
print "// PointInPrimaryASU: function for determining whether a given point lies within",
print "the primary"
print "//                    asymmetric unit as defined by the space group's geometry."
print "//"
print "// Arguments:"
print "//   [x,y,z]:    the fractional coordinates of the point"
print "//" + ("-" * 93)
print "bool Action_XtalSymm::PointInPrimaryASU(double x, double y, double z)"
print "{"
print "  const double two3 = 0.66666666666667;"
print "  const double half = 0.5;"
print "  const double thr8 = 0.375;"
print "  const double third = 0.33333333333333;"
print "  const double fourth = 0.25;"
print "  const double sixth = 0.16666666666667;"
print "  const double eighth = 0.125;"
print "  const double twelfth = 0.83333333333333;"
print "  bool result = true;"
print "  switch(sgID) {"
for i, item in enumerate(grpcomps):
  print "    case %d:" % i
  lsize = 0
  noreset = True;
  for j, comparison in enumerate(item):
    lhs = comparison[0]
    lhs = lhs.replace("2/3",  "twothr")
    lhs = lhs.replace("1/2",  "half")
    lhs = lhs.replace("3/8",  "threig")
    lhs = lhs.replace("1/3",  "third")
    lhs = lhs.replace("1/4",  "fourth")
    lhs = lhs.replace("1/6",  "sixth")
    lhs = lhs.replace("1/8",  "eighth")
    lhs = lhs.replace("1/12", "twelfth")
    lhs = lhs.replace("0", "0.0")
    lhs = lhs.replace("1", "1.0")
    lhs = lhs.replace("2", "2.0")
    lhs = lhs.replace("3", "3.0")
    lhs = lhs.replace(",", ", ")
    lhs = lhs.replace(".0x", ".0*x")
    lhs = lhs.replace(".0y", ".0*y")
    lhs = lhs.replace(".0z", ".0*z")
    lhs = lhs.replace("min", "dmin");
    lhs = lhs.replace("max", "dmax");
    rhs = comparison[1]
    rhs = rhs.replace("2/3",  "twothr")
    rhs = rhs.replace("1/2",  "half")
    rhs = rhs.replace("3/8",  "threig")
    rhs = rhs.replace("1/3",  "third")
    rhs = rhs.replace("1/4",  "fourth")
    rhs = rhs.replace("1/6",  "sixth")
    rhs = rhs.replace("1/8",  "eighth")
    rhs = rhs.replace("1/12", "twelfth")
    rhs = rhs.replace("0", "0.0")
    rhs = rhs.replace("1", "1.0")
    rhs = rhs.replace("2", "2.0")
    rhs = rhs.replace("3", "3.0")
    rhs = rhs.replace(",", ", ")
    rhs = rhs.replace(".0x", ".0*x")
    rhs = rhs.replace(".0y", ".0*y")
    rhs = rhs.replace(".0z", ".0*z")
    rhs = rhs.replace("min", "dmin");
    rhs = rhs.replace("max", "dmax");
    if (j == 0):
      contrib = "      if (%s > %s" % (lhs, rhs)
      if (j < len(item) - 1):
        contrib += " ||"
      else:
        contrib += ")"
      lsize = len(contrib) + 1
      print contrib,
    else:
      contrib = "%s > %s" % (lhs, rhs)
      if (j < len(item) - 1):
        contrib += " ||"
      else:
        contrib += ")"
      if (lsize + len(contrib) < 92):
        lsize += len(contrib) + 1
        print contrib,
      else:
        print
        noreset = False;
        contrib = "          " + contrib
        lsize = len(contrib) + 1
        print contrib,
  contrib = "return false;"
  if (lsize + len(contrib) < 96 and noreset):
    print contrib
  else:
    print "{"
    print "        return false;"
    print "      }"
print "  }"
print "  return result;"
print "}"
