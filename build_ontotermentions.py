import os
import argparse
import sys

if __name__ == '__main__':
    print("File one executed when ran directly")
    parser=argparse.ArgumentParser()
    parser.add_argument('--path', '-i',help='Name of the path')
    
    args=parser.parse_args()
    path = args.path
    
    
    if args is None :
        parser.print_help()
        sys.exit('Not enough arguments. Expected at least -i "Path to ontology folder" ')
    else:
        print(path)
    
    if path is None: 
        parser.print_help()
        sys.exit('Not enough arguments. Expected at least -i "Path to ontology folder" ')
    else:
        print("got it: ", path)
        #todo: get ontopath from current directory
        ontopath = os.path.dirname(os.path.realpath(__file__))
        print("ontopath is: ", ontopath)
        # ontopath = "/home/tom/Documents/PROGRAMMING/Python/Flask/onto-text-tag"

        addictopath = path + "/scripts/"
        # os.chdir(addictopath)
        # exec(open("./1getPubMedAbstractsAuthors.py").read())
        os.chdir(ontopath)
        exec(open("./ontotagtext.py").read())
        os.chdir(addictopath)
        exec(open("./3applyTextTaggingToAbstracts.py").read())

    # os.chdir("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/")
    # exec(open("./4buldDataFrame.py").read())
    # os.chdir("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/")
    # exec(open("./5otherGraphs.py").read())

    # os.chdir("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/")
    # exec(open("./6buildOntoBio.py").read())
    # print("df is: ", df)
    # df.to_csv("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/ontotermentins.csv", sep = ",", encoding= "utf-8")
