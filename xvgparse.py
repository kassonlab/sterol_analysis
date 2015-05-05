#!/usr/bin/python

#Break out xvg parser
#by Peter Kasson

def xvg_parse( infile, columns, keycol=-1, skipcol=[], duplicates=[], ascending=1 ):
  #utility routine to deal with xvg output files (Gromacs tools)
  """ From FAHAnalysis Tools
      code to parse an xvg file
      Passing None to columns indicates include all columns not in skipcol
  """
  f = open( infile )
  valuelist = {} # keys from keycol

  numLines = 0
  for line in f:
    line = line.strip()
    values = []

    #skip comments and blank lines and convert to numbers
    if line and ( line[0] not in ['@', '#'] ):
      numLines += 1
      words = line.split()
      # get values
      try:
        if columns is None:
          for col in range(len(words)):
            if col not in skipcol:
              values.append( float( words[col] ) )
        else:
          for col in columns:
            values.append( float( words[col] ) )
        if keycol != -1:
          keyval = float(words[keycol])
        else:
          keyval = numLines
        # if ascending and find duplicate use most recent value
        if ascending == 1:
          if valuelist.has_key(keyval):
            print "-- duplicate snapshot in xvg, overwriting with new value at line " + str(numLines) 
          valuelist[keyval] = values
        # if descending and find duplicate then don'e overwrite
        elif ascending == -1:
          if valuelist.has_key(keyval):
            print "-- duplicate snapshot in xvg, skipping line " + str(numLines)
          else:
            valuelist[keyval] = values
        else:
          valuelist[keyval] = values
      except Exception:
        # print "-- Corrupted xvg file"
        # return []
        print "could not parse xvg line:", words
        numLines -= 1
        continue

  # create list to return
  # should be ordered by keys, but not include keys themselves
  returnlist = []
  keylist = valuelist.keys()
  keylist.sort()
  for mykey in keylist:
    returnlist.append(valuelist[mykey])

  return returnlist

    
