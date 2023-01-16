import re
import codecs
import os
import argparse

# Function to get items from the File
def getInputs(line, inputs):
    flag = 0
    # Take the inputs from the line
    for i in range(3, len(line.split(" "))-2, 2):
        if inputs[len(inputs)-1] == '' or flag:
            inputs.append(line.split(" ")[i+1].split(",")[0].strip())
            flag = 1
        else:
            inputs.append(line.split(" ")[i].split(",")[0].strip())
    # Check if the number of threads is defined or not
    if 'OMP_NUM_THREADS' in locals():
        command = "OMP_NUM_THREADS=" + n_threads + " ./test_lapack.x " + api_name + " "
    else:
        command = "./test_lapack.x " + api_name + " "
    # Append the inputs
    for param in inputs:
        if param == '':
            command = command + "  "
        else:
            command = command + param + " "
    # Append the repeats at the end
    command = command + str(args.nrepeats)
    commands.add(command)


# Global Variable
commands = set()

# Check if Number of threads is defined
if os.environ.get('OMP_NUM_THREADS') is not None:
    n_threads = os.environ['OMP_NUM_THREADS']

# Take filename,api name and repeats from the arguments
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True)
parser.add_argument("--nrepeats", type=int, default=1)
parser.add_argument("--apiname", type=str)
args = parser.parse_args()

# Read the log file and fetch the api and the inputs
with codecs.open(args.filename, 'r', encoding='utf-8', errors='ignore') as file_log:
    for line in file_log:
        # Skip the line if it has tabs/space
        if (re.search("^(\t+|\s+)", line) == None):
            inputs = []
            # Split the api name from the line and the precision form the api name
            api_name = line.split(" ")[0].split(":")[2]
            inputs.append(api_name[0])
            api_name = api_name[1:]
            # Check if the apiName is passed in the arguments
            if args.apiname is not None:
                if args.apiname == api_name:
                    getInputs(line,inputs)
            else:
                getInputs(line,inputs)
file_log.close()

# Execute all the apis
os.chdir('..')
for command in commands:
    os.system(command)
