from sys import *

def main():
    
    # Open files
    inf = open("../tables/%s_%s_%s_%s.txt" % (argv[1], argv[2], argv[3], argv[4]), "r")
    out = open("%s_%s_%s_%s.csv" % (argv[1], argv[2], argv[3], argv[4]), "w")
    out.write("N,t,E,D,D1,E2,D2,E3,D3 \n")

    # Process file line by line
    for line in inf:
        
        line = line.strip()
        line = line.strip("\\")
        
        line = line.split("&")
        line = list(filter(None, line))
        
        for value in line:
            out.write(value)
            out.write(",")
        out.write("\n")

    # Close files
    inf.close()
    out.close()

main()
