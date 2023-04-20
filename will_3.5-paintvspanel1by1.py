import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-idfile",
                    help="file containing all IDs to process; only those in pop 'UKBB' will be processed in current implementation",
                    required=True)
parser.add_argument("-cp_panel_scripts",
                    help="the location of the panel scripts which will be copied to each temp directory",
                    required=True)
args = parser.parse_args()
idfile = pd.read_csv(args.idfile, sep=" ", header=None)
idfile = idfile[idfile[1]=='UKBB']
idlist = idfile[0].tolist()

myfile = open('paintvspanel1by1_commands.txt', 'w')
for n,id in enumerate(idlist):
    myfile.write("bash paintsample1by1.sh " + str(id) + " " + str(n + 1) + " " + str((n+2)*2) + " " + args.cp_panel_scripts + "\n")
myfile.close()

