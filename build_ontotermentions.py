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
        #Get ontopath from current directory
        ontopath = os.path.dirname(os.path.realpath(__file__))
        print("ontopath is: ", ontopath)
        # ontopath = "/home/tom/Documents/PROGRAMMING/Python/Flask/onto-text-tag"

        addictopath = path + "/scripts/"
        
        # #delete all the db files in the DB folder:
        os.chdir(addictopath + "DB/")
        for the_file in os.listdir(addictopath + "DB/"):
            file_path = os.path.join(addictopath + "DB/", the_file)
            try:
                if os.path.isfile(file_path):
                    if file_path.endswith('.bak') or file_path.endswith('.dat') or file_path.endswith('.dir'): 
                        os.unlink(file_path)
            except Exception as e:
                print(e)

        os.chdir(addictopath)
        exec(open("./1getPubMedAbstractsAuthors.py").read())
        os.chdir(ontopath)
        exec(open("./ontotagtext.py").read())
        os.chdir(addictopath)
        exec(open("./3applyTextTaggingToAbstracts.py").read())

    #delete all db files in the static folder:
        os.chdir(ontopath + "/static/")
        for the_file in os.listdir(ontopath + "/static"):
            file_path = os.path.join(ontopath + "/static/", the_file)
            try:
                if os.path.isfile(file_path):
                    if file_path.endswith('.bak') or file_path.endswith('.dat') or file_path.endswith('.dir'): 
                        os.unlink(file_path)
            except Exception as e:
                print(e)

    #move all files in the DB folder to the ontology folder:
    os.chdir(addictopath + "DB/")
    for the_file in os.listdir(addictopath + "DB/"):
        file_path = os.path.join(addictopath + "DB/", the_file)
        try:
            if os.path.isfile(file_path):
                if file_path.endswith('.bak') or file_path.endswith('.dat') or file_path.endswith('.dir'): 
                    os.rename(file_path, ontopath + "/static/" + the_file)
        except Exception as e:
            print(e)
       



    #OLD CODE:
    # os.chdir("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/")
    # exec(open("./4buldDataFrame.py").read())
    # os.chdir("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/")
    # exec(open("./5otherGraphs.py").read())

    # os.chdir("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/")
    # exec(open("./6buildOntoBio.py").read())
    # print("df is: ", df)
    # df.to_csv("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/ontotermentins.csv", sep = ",", encoding= "utf-8")
