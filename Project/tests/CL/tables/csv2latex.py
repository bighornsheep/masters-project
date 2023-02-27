from sys import *

def main():
    
    # Open files
    inf = open("../csv_16_80_256/%s_%s_%s_%s.csv" % (argv[1], argv[2], argv[3], argv[4]), "r")
    out = open("%s_%s_%s_%s.txt"        % (argv[1], argv[2], argv[3], argv[4]), "w")

    lines = inf.readlines()[1:]

    # Process file line by line
    for line in lines:
        line = line.strip("\n")
        line = line.split(",")
        line = list(filter(None, line))

        
        for i in range (len(line)):
            out.write (line[i])
            if i != len(line) - 1:
                out.write (" & ")
        out. write (" \\\\")
        out. write (" \n")
        
        
    # Close files
    inf.close()
    out.close()

main()
